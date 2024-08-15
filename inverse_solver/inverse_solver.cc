#include "inverse_solver.h"
#include "opensn/modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "opensn/modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_linear_solver.h"
#include "opensn/framework/mesh/mesh_continuum/mesh_continuum.h"
#include "opensn/framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "opensn/framework/field_functions/field_function_grid_based.h"
#include "opensn/framework/materials/multi_group_xs/multi_group_xs.h"
#include "opensn/framework/logging/log.h"
#include "opensn/framework/runtime.h"
#include "opensn/framework/object_factory.h"
#include <numeric>
#include <iomanip>

namespace opensn
{

std::string
TaoConvergedReasonString(TaoConvergedReason reason)
{
  std::stringstream ss;
  switch (reason)
  {
    // Converged
    case TAO_CONVERGED_GATOL:
      ss << "TAO_CONVERGED_GATOL";
      break;
    case TAO_CONVERGED_GRTOL:
      ss << "TAO_CONVERGED_GRTOL";
      break;
    case TAO_CONVERGED_GTTOL:
      ss << "TAO_CONVERGED_GTTOL";
      break;
    case TAO_CONVERGED_STEPTOL:
      ss << "TAO_CONVERGED_STEPTOL";
      break;
    case TAO_CONVERGED_MINF:
      ss << "TAO_CONVERGED_MINF";
      break;
    case TAO_CONVERGED_USER:
      ss << "TAO_CONVERGED_USER";
      break;

    // Diverged
    case TAO_DIVERGED_MAXITS:
      ss << "TAO_DIVERGED_MAXITS";
      break;
    case TAO_DIVERGED_NAN:
      ss << "TAO_DIVERGED_NAN";
      break;
    case TAO_DIVERGED_MAXFCN:
      ss << "TAO_DIVERGED_MAXFCN";
      break;
    case TAO_DIVERGED_LS_FAILURE:
      ss << "TAO_DIVERGED_LS_FAILURE";
      break;
    case TAO_DIVERGED_TR_REDUCTION:
      ss << "TAO_DIVERGED_TR_REDUCTION";
      break;
    case TAO_DIVERGED_USER:
      ss << "TAO_DIVERGED_USER";
      break;

    // Unknown
    default:
      ss << "Unknown TaoConvergedReason";
      break;
  }
  return ss.str();
}

PetscErrorCode
__EvaluateObjectiveAndGradient(Tao, Vec x, PetscReal* f, Vec df, void* ctx)
{
  auto app_ctx = static_cast<InverseSolver*>(ctx);
  return app_ctx->EvaluateObjectiveAndGradient(x, f, df);
}

PetscErrorCode
__Monitor(Tao, void* ctx)
{
  auto app_ctx = static_cast<InverseSolver*>(ctx);
  return app_ctx->Monitor();
}

OpenSnRegisterObjectInNamespace(lbs, InverseSolver);

InputParameters
InverseSolver::GetInputParameters()
{
  auto params = Solver::GetInputParameters();

  params.SetGeneralDescription("An inverse solver for density reconstruction.");
  params.ChangeExistingParamToOptional("name", "LBSInverseDataBlock");

  params.AddRequiredParameter<size_t>("lbs_solver_handle", "A handle to an LBS solver.");
  params.AddRequiredParameterArray("forward_bcs",
                                   "A parameter block for forward boundary conditions.");

  params.AddRequiredParameterArray(
    "detector_boundaries",
    "An array of detector boundary names. These boundaries represent those included in the "
    "minimization functional.The available boundary names are \"xmin\", \"xmax\", \"ymin\", "
    "\"ymax\", \"zmin\", and \"zmax\".");

  params.AddOptionalParameter(
    "material_ids", std::vector<int>(), "The material IDs with unknown densities.");
  params.AddOptionalParameter(
    "initial_guess", std::vector<double>(), "An initial guess for material-based densities.");
  params.AddOptionalParameter(
    "initial_guess_file", "", "A file with the initial guess for cell-based densities.");
  params.AddOptionalParameter(
    "cellwise", false, "A flag to solve the inverse problem for cell-wise densities.");

  params.AddOptionalParameter("max_its", 20, "The maximum number of iterations.");
  params.AddOptionalParameter("tol", 1.0e-8, "The convergence tolerance.");
  params.AddOptionalParameter("alpha", 1.0, "The initial step size.");
  params.AddOptionalParameter(
    "line_search", false, "A flag for using a line search algorithm for step length selection.");
  params.AddOptionalParameter("max_ls_its", 20, "The maximum number of line search iterations.");
  params.AddOptionalParameter("use_tao", false, "A flag to use TAO from PETSc.");

  params.ConstrainParameterRange("max_its", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("tol", AllowableRangeLowLimit::New(1.0e-16));
  params.ConstrainParameterRange("alpha", AllowableRangeLowHighLimit::New(1.0e-3, 1.0e6));
  params.ConstrainParameterRange("max_ls_its", AllowableRangeLowLimit::New(1));

  return params;
}

InverseSolver::InverseSolver(const InputParameters& params)
  : Solver(params),
    solver_(GetStackItem<DiscreteOrdinatesSolver>(
      object_stack, params.GetParamValue<size_t>("lbs_solver_handle"))),
    forward_bcs_(params.GetParam("forward_bcs")),
    detector_bndrys_(params.GetParamVectorValue<std::string>("detector_boundaries")),
    num_transport_solves_(0),
    max_its_(params.GetParamValue<unsigned>("max_its")),
    tol_(params.GetParamValue<double>("tol")),
    alpha_(params.GetParamValue<double>("alpha")),
    line_search_(params.GetParamValue<bool>("line_search")),
    max_ls_its_(params.GetParamValue<unsigned>("max_ls_its")),
    use_tao_(params.GetParamValue<bool>("use_tao")),
    cellwise_(params.GetParamValue<bool>("cellwise"))
{
  const auto user_params = params.ParametersAtAssignment();

  // Check boundary condition options
  forward_bcs_.RequireBlockTypeIs(ParameterBlockType::ARRAY);

  // Initial guess set from material IDs
  if (user_params.Has("material_ids"))
  {
    user_params.RequireParameter("material_ids");
    user_params.RequireParameter("initial_guess");

    material_ids_ = user_params.GetParamVectorValue<int>("material_ids");
    const auto initial_guess = user_params.GetParamVectorValue<double>("initial_guess");
    OpenSnLogicalErrorIf(initial_guess.size() != material_ids_.size(),
                         "The initial guess must have the same number of entries as the "
                         "specified material IDs.");

    // Define the initial guess vector. If cell-wise, map the material-based densities to
    // a vector of densities for cells containing the given materials.
    std::vector<double> rho0(initial_guess);
    if (cellwise_)
    {
      rho0.clear();
      for (const auto& cell : GetCurrentMesh()->local_cells)
      {
        const auto& it = std::find(material_ids_.begin(), material_ids_.end(), cell.material_id_);
        if (it != material_ids_.end())
        {
          cell_local_ids_.push_back(cell.local_id_);
          rho0.push_back(initial_guess[it - material_ids_.begin()]);
        }
      }
    }

    // Set the initial guess
    VecCreate(opensn::mpi_comm, &solution_);
    VecSetType(solution_, VECSTANDARD);
    VecSetSizes(solution_, rho0.size(), PETSC_DECIDE);

    double* val;
    VecGetArray(solution_, &val);
    for (PetscInt i = 0; i < rho0.size(); ++i)
      val[i] = rho0[i];
    VecRestoreArray(solution_, &val);

    // Create the gradient vector
    VecDuplicate(solution_, &gradient_);
    VecSet(gradient_, 0.0);
  }
  else
  {
    throw std::invalid_argument("Only material-based initialization is implemented.");
  }

  if (use_tao_)
  {
    TaoCreate(opensn::mpi_comm, &tao_);
    TaoSetType(tao_, TAOBLMVM);
    TaoSetSolution(tao_, solution_);
    TaoSetApplicationContext(tao_, this);
    TaoSetObjectiveAndGradient(tao_, gradient_, __EvaluateObjectiveAndGradient, this);
    TaoSetMonitor(tao_, __Monitor, this, NULL);

    TaoSetMaximumIterations(tao_, max_its_);
    TaoSetTolerances(tao_, tol_, tol_, tol_);
    TaoSetFromOptions(tao_);
  }
}

InverseSolver::~InverseSolver()
{
  VecDestroy(&solution_);
  VecDestroy(&gradient_);
  if (use_tao_)
    TaoDestroy(&tao_);
}

void
InverseSolver::Initialize()
{
  // Check for invalid material IDs
  if (not material_ids_.empty())
    for (const auto& matid : material_ids_)
      OpenSnLogicalErrorIf(solver_.GetMatID2XSMap().count(matid) == 0,
                           "Material ID " + std::to_string(matid) + " in the density map " +
                             "does not correspond to a material ID in the simulation.");

  // Define the detector boundary IDs from names
  detector_bids_.reserve(detector_bndrys_.size());
  for (const auto& bname : detector_bndrys_)
    detector_bids_.push_back(solver_.supported_boundary_names.at(bname));

  // Generate the measured signal to reconstruct.
  // NOTE The densities in the solver should be set to the nominal solution at this point.
  solver_.SetOptions(GetForwardOptions());
  ExecuteSteadyState();
  ref_leakage_ = ComputeDetectorSignal();

  log.Log() << "\n********** Reference Signal **********\n"
               "Boundary         Value\n"
               "--------         -----\n";
  for (int i = 0; i < ref_leakage_.size(); ++i)
    log.Log() << std::left << std::setw(8) << detector_bndrys_[i] << "         " << std::left
              << std::setw(8) << std::setprecision(8) << ref_leakage_[i];
  log.Log() << "\n";

  // Set densities to the initial guess
  SetDensities(solution_, true);

  // Reset xs to microscopic
  for (const auto& [matid, xs] : solver_.GetMatID2XSMap())
    std::dynamic_pointer_cast<MultiGroupXS>(xs)->SetScalingFactor(1.0);
}

void
InverseSolver::Execute()
{
  if (use_tao_)
    ExecuteTao();
  else
    ExecuteSimple();

  log.Log() << "\n********** Inverse Solver Summary **********";
  log.Log() << "Solution:            " << VecString(solution_);
  log.Log() << "# Transport Solves:  " << num_transport_solves_;
  log.Log() << "\n";

  SetDensities(solution_);
  ExportDensity();
}

void
InverseSolver::ExecuteTao()
{
  TaoSolve(tao_);
  log.Log() << "\n********** Tao Solver Summary **********";
  TaoView(tao_, PETSC_VIEWER_STDOUT_WORLD);
  log.Log() << "\n";
}

void
InverseSolver::ExecuteSimple()
{
  // Bookkeeping
  VecDuplicate(solution_, &gradient_);
  VecSet(gradient_, 0.0);

  // Start iterations
  double f, alpha;
  for (int nit = 0; nit < max_its_; ++nit)
  {
    EvaluateObjectiveAndGradient(solution_, &f, gradient_);
    alpha = line_search_ ? BackTrackingLineSearch(solution_, alpha_, gradient_, f) : alpha_;
    VecAXPY(solution_, -alpha, gradient_);

    if (nit == 0)
      log.Log() << "********** Iteration Summary **********";

    std::stringstream ss;
    ss.precision(6);
    ss << "Iteration " << std::right << std::setw(4) << nit << "  ";
    ss << "Obj  " << std::right << std::setw(12) << f << "  ";
    ss << "Alpha  " << std::right << std::setw(12) << alpha << "  ";
    ss << "Solution  " << VecString(solution_) << "  ";
    ss << "Gradient  " << VecString(gradient_) << (f < tol_ ? "    CONVERGED" : "");
    log.Log() << ss.str();

    if (f < tol_)
      break;
  }
  log.Log() << "\n********** " << num_transport_solves_ << " transport solves";
  log.Log() << "\n";
}

PetscErrorCode
InverseSolver::EvaluateObjective(Vec x, PetscReal* f)
{
  // Solve the forward problem
  SetDensities(x);
  solver_.SetOptions(GetForwardOptions());
  ExecuteSteadyState();

  // Evaluate the objective function
  *f = 0.0;
  const auto leakage = ComputeDetectorSignal();
  for (int i = 0; i < leakage.size(); ++i)
    *f += 0.5 * std::pow(leakage[i] - ref_leakage_[i], 2.0);
  return 0;
}

PetscErrorCode
InverseSolver::EvaluateObjectiveAndGradient(Vec x, PetscReal* f, Vec df)
{
  // Solve the forward problem
  SetDensities(x);
  solver_.SetOptions(GetForwardOptions());
  ExecuteSteadyState();

  // Evaluate the objective function
  *f = 0.0;
  const auto leakage = ComputeDetectorSignal();
  for (int i = 0; i < leakage.size(); ++i)
    *f += 0.5 * std::pow(leakage[i] - ref_leakage_[i], 2.0);

  // Get data for the inner product
  const auto phi_fwd = solver_.PhiNewLocal();
  const auto psi_fwd = solver_.PsiNewLocal();

  // Solve the adjoint problem
  solver_.SetOptions(GetAdjointOptions(leakage));
  ExecuteSteadyState();
  const auto psi_adj = solver_.PsiNewLocal();

  // Compute the gradient
  ComputeInnerProduct(phi_fwd, psi_fwd, psi_adj, df);
  return 0;
}

PetscErrorCode
InverseSolver::Monitor() const
{
  Vec x, g;
  PetscInt its, n;
  PetscReal f, gnorm, cnorm, xdiff;
  TaoConvergedReason reason;

  TaoGetSolution(tao_, &x);
  VecGetLocalSize(x, &n);
  TaoGetGradient(tao_, &g, NULL, NULL);
  TaoGetSolutionStatus(tao_, &its, &f, &gnorm, &cnorm, &xdiff, &reason);

  std::stringstream ss;
  ss << its << ",  Function Value: " << f << ",  Gradient Norm: " << gnorm
     << ",  Step Length: " << xdiff;
  if (n <= 10)
    ss << ",  Solution: " << VecString(x) << ",  Gradient: " << VecString(g);
  log.Log() << ss.str();

  return 0;
}

std::vector<double>
InverseSolver::ComputeDetectorSignal() const
{
  std::vector<double> leakage;
  leakage.reserve(detector_bids_.size());
  const auto leakage_map = solver_.ComputeLeakage(detector_bids_);
  for (const auto& [bid, vals] : leakage_map)
    leakage.push_back(std::accumulate(vals.begin(), vals.end(), 0.0));
  return leakage;
}

void
InverseSolver::ComputeInnerProduct(const std::vector<double>& phi_fwd,
                                   const std::vector<std::vector<double>>& psi_fwd,
                                   const std::vector<std::vector<double>>& psi_adj,
                                   Vec df) const
{
  // Zero out the df vector
  VecSet(df, 0.0);

  // Set to forward mode
  solver_.SetOptions(GetForwardOptions());

  // Get solver data
  const auto& grid = solver_.Grid();
  const auto& discretization = solver_.SpatialDiscretization();
  const auto& unit_cell_matrices = solver_.GetUnitCellMatrices();
  const auto& transport_views = solver_.GetCellTransportViews();
  const auto num_moments = solver_.NumMoments();

  // Loop over groupsets
  for (const auto& groupset : solver_.Groupsets())
  {
    const auto gs = groupset.id_;
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature_;
    const auto& m2d = quadrature->GetMomentToDiscreteOperator();
    const auto& moment_map = quadrature->GetMomentToHarmonicsIndexMap();
    const auto& weights = quadrature->weights_;
    const auto num_gs_dirs = quadrature->omegas_.size();
    const auto num_gs_groups = groupset.groups_.size();
    const auto first_gs_group = groupset.groups_.front().id_;

    // Loop over cells
    for (const auto& cell : grid.local_cells)
    {
      // Default to cell-wise index
      const auto idx = GetIndex(cell);
      if (idx < 0)
        continue;

      // Get cell data
      const auto& cell_mapping = discretization.GetCellMapping(cell);
      const auto& transport_view = transport_views[cell.local_id_];
      const auto& mass_matrix = unit_cell_matrices[cell.local_id_].intV_shapeI_shapeJ;
      const auto num_cell_nodes = cell_mapping.NumNodes();
      const auto& xs = transport_views[cell.local_id_].XS();
      const auto& rho = cellwise_ ? solver_.DensitiesLocal()[cell.local_id_] : xs.ScalingFactor();

      // Get xs data
      const auto& sigma_t = xs.SigmaTotal();
      const auto& S = xs.TransferMatrices();
      const auto& F = xs.ProductionMatrix();

      // Loop over cell nodes
      double val = 0.0;
      for (int i = 0; i < num_cell_nodes; ++i)
      {
        for (int j = 0; j < num_cell_nodes; ++j)
        {
          const auto m_ij = mass_matrix[i][j];

          // Loop over groupset groups
          for (int gsg = 0; gsg < num_gs_groups; ++gsg)
          {
            const auto g = first_gs_group + gsg;
            const auto& sig_t = sigma_t[g] / rho;

            // Precompute source moments
            std::vector<double> q_mom(num_moments, 0.0);
            for (int m = 0; m < num_moments; ++m)
            {
              const auto ell = moment_map[m].ell;
              const auto uk_map = transport_view.MapDOF(j, m, 0);

              // Scattering
              double s = 0.0;
              if (ell < S.size())
                for (const auto& [_, gp, sig_ell] : S[ell].Row(g))
                  s += sig_ell * phi_fwd[uk_map + gp];

              // Fission
              double f = 0.0;
              unsigned gp = 0;
              if (xs.IsFissionable() and ell == 0)
                for (const auto& sig_f : F[g])
                  f += sig_f * phi_fwd[uk_map + gp++];

              // Accumumlate source moment
              q_mom[m] = (s + f) / rho;
            } // for moment

            // Angular integration
            for (int d = 0; d < num_gs_dirs; ++d)
            {
              const auto w = weights[d] * m_ij;
              const auto dof_adj = discretization.MapDOF(cell, i, uk_man, d, gsg);
              const auto dof_fwd = discretization.MapDOF(cell, j, uk_man, d, gsg);

              // Apply moment-to-discrete operator to source moments
              double q_dir = 0.0;
              for (int m = 0; m < num_moments; ++m)
                q_dir += m2d[m][d] * q_mom[m];

              // Inner product contribution
              val += w * psi_adj[gs][dof_adj] * (sig_t * psi_fwd[gs][dof_fwd] - q_dir);
            } // for direction
          }   // for group g
        }     // for node j
      }       // for node i

      VecSetValue(df, idx, -val, ADD_VALUES);

    } // for cell
  }   // for groupset

  VecAssemblyBegin(df);
  VecAssemblyEnd(df);
}

double
InverseSolver::BackTrackingLineSearch(Vec x, double alpha0, Vec df, double f0)
{
  Vec xp;
  VecDuplicate(x, &xp);

  // Define the sufficient decrease critereia
  double m = 0.0;
  VecNorm(df, NORM_2, &m);
  m = -m * m;
  double t = -1.0e-4 * m;

  // Start the line search
  double alpha(alpha0), f(f0);
  for (int i = 0; i < max_ls_its_; ++i)
  {
    VecCopy(x, xp);
    VecAXPY(xp, -alpha, df);
    EvaluateObjective(xp, &f);
    if (f0 - f >= alpha * t)
      return alpha;
    alpha *= 0.5;
  }

  VecDestroy(&xp);
  return alpha;
}

void
InverseSolver::SetDensities(Vec x, const bool init)
{
  int64_t n;
  VecGetLocalSize(x, &n);

  if (cellwise_ and n != cell_local_ids_.size())
    throw std::logic_error("Vector size mismatch.");
  if (not cellwise_ and n != material_ids_.size())
    throw std::logic_error("Vector size mismatch.");

  const double* val;
  VecGetArrayRead(x, &val);
  auto& densities = solver_.DensitiesLocal();
  for (const auto& cell : solver_.Grid().local_cells)
  {
    const auto idx = GetIndex(cell);
    if (idx >= 0)
      densities[cell.local_id_] = val[idx];
    else if (init)
      densities[cell.local_id_] = solver_.GetMatID2XSMap().at(cell.material_id_)->ScalingFactor();
  }
  VecRestoreArrayRead(x, &val);
}

int
InverseSolver::GetIndex(const Cell& cell) const
{
  if (not cellwise_)
  {
    const auto it = std::find(material_ids_.begin(), material_ids_.end(), cell.material_id_);
    return it != material_ids_.end() ? it - material_ids_.begin() : -1;
  }
  else
  {
    const auto it = std::find(cell_local_ids_.begin(), cell_local_ids_.end(), cell.local_id_);
    return it != cell_local_ids_.end() ? it - cell_local_ids_.begin() : -1;
  }
}

InputParameters
InverseSolver::GetForwardOptions() const
{
  ParameterBlock params;
  params.AddParameter("adjoint", false);
  params.AddParameter("clear_boundary_conditions", true);

  ParameterBlock bcs("boundary_conditions");
  for (int b = 0; b < forward_bcs_.NumParameters(); ++b)
    bcs.AddParameter(forward_bcs_.GetParam(b));
  bcs.ChangeToArray();
  params.AddParameter(bcs);

  auto lbs_params = LBSSolver::OptionsBlock();
  lbs_params.AssignParameters(params);
  return lbs_params;
}

InputParameters
InverseSolver::GetAdjointOptions(const std::vector<PetscReal>& leakage) const
{
  ParameterBlock params;
  params.AddParameter("adjoint", true);
  params.AddParameter("clear_boundary_conditions", true);

  ParameterBlock bcs("boundary_conditions");
  for (int i = 0; i < detector_bndrys_.size(); ++i)
  {
    const auto mismatch = leakage[i] - ref_leakage_[i];
    std::vector<double> vals(solver_.NumGroups(), mismatch);

    ParameterBlock bc(std::to_string(i + 1));
    bc.AddParameter("name", detector_bndrys_[i]);
    bc.AddParameter("type", "isotropic");
    ParameterBlock group_strength("group_strength", vals);
    bc.AddParameter(group_strength);
    bcs.AddParameter(bc);
  }
  bcs.ChangeToArray();
  params.AddParameter(bcs);

  auto lbs_params = LBSSolver::OptionsBlock();
  lbs_params.AssignParameters(params);
  return lbs_params;
}

void
InverseSolver::ExecuteSteadyState()
{
  ++num_transport_solves_;

  auto& ags_solver = *solver_.GetPrimaryAGSSolver();
  ags_solver.Setup();
  ags_solver.Solve();

  if (solver_.Options().use_precursors)
    solver_.ComputePrecursors();

  if (solver_.Options().adjoint)
    solver_.ReorientAdjointSolution();

  solver_.UpdateFieldFunctions();
}

void
InverseSolver::ExportDensity() const
{
  const auto& grid = solver_.Grid();
  const auto& discretization = solver_.SpatialDiscretization();
  const auto& cs_type = discretization.GetCoordinateSystemType();
  const auto& density = solver_.DensitiesLocal();

  auto fv = std::dynamic_pointer_cast<SpatialDiscretization>(FiniteVolume::New(grid, cs_type));
  std::vector<std::shared_ptr<const FieldFunctionGridBased>> ff_list;
  ff_list.emplace_back(std::make_shared<const FieldFunctionGridBased>(
    "density", fv, Unknown(UnknownType::SCALAR), density));
  FieldFunctionGridBased::ExportMultipleToVTK("density", ff_list);
}

std::string
InverseSolver::VecString(Vec x)
{
  int64_t n;
  const double* ptr;

  std::stringstream ss;
  ss.precision(6);

  VecGetLocalSize(x, &n);
  VecGetArrayRead(x, &ptr);
  ss << "[";
  for (int64_t i = 0; i < n - 1; ++i)
    ss << ptr[i] << "  ";
  ss << ptr[n - 1] << "]";
  VecRestoreArrayRead(x, &ptr);

  return ss.str();
}

} // namespace opensn
