#include "lbs_inverse_solver.h"
#include "opensn/modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_linear_solver.h"
#include "opensn/framework/mesh/mesh_continuum/mesh_continuum.h"
#include "opensn/framework/physics/physics_material/multi_group_xs/single_state_mgxs.h"
#include "opensn/framework/logging/log.h"
#include "opensn/framework/runtime.h"
#include "opensn/framework/object_factory.h"
#include <numeric>
#include <iomanip>

namespace opensn::lbs
{

OpenSnRegisterObjectInNamespace(lbs, InverseSolver);

InputParameters
InverseSolver::GetInputParameters()
{
  auto params = Solver::GetInputParameters();

  params.SetGeneralDescription("An inverse solver for material density reconstruction.");

  params.ChangeExistingParamToOptional("name", "LBSInverseDataBlock");

  params.AddRequiredParameter<size_t>("lbs_solver_handle",
                                      "A handle to the solver used in the inverse problem.");

  params.AddRequiredParameterArray(
    "boundary_conditions", "A parameter block used to set the forward boundary conditions.");

  params.AddRequiredParameterArray(
    "detector_boundaries",
    "An array of detector boundary names. These boundaries represent those included in the "
    "minimization functional.The available boundary names are \"xmin\", \"xmax\", \"ymin\", "
    "\"ymax\", \"zmin\", and \"zmax\".");

  params.AddRequiredParameterBlock(
    "initial_guess",
    "A table containing an array named \"material_ids\" containing the "
    "material IDs to be reconstructed and an array named \"values\" "
    "containing the initial guesses per material ID.");

  params.AddOptionalParameter("alpha", 1.0, "The initial step size.");

  params.AddOptionalParameter("max_iterations", 10, "The maximum number of iterations to take.");

  params.AddOptionalParameter(
    "tolerance", 1.0e-6, "The convergence tolerance for the density reconstruction.");

  params.AddOptionalParameter("line_search", true, "A flag for using a line search algorithm.");

  return params;
}

InverseSolver::InverseSolver(const InputParameters& params)
  : Solver(params),
    solver_(GetStackItem<DiscreteOrdinatesSolver>(
      object_stack, params.GetParamValue<size_t>("lbs_solver_handle"))),
    bc_options_(params.GetParam("boundary_conditions")),
    detector_boundaries_(params.GetParamVectorValue<std::string>("detector_boundaries")),
    alpha_max_(params.GetParamValue<double>("alpha")),
    max_iterations_(params.GetParamValue<unsigned int>("max_iterations")),
    tolerance_(params.GetParamValue<double>("tolerance")),
    line_search_(params.GetParamValue<bool>("line_search"))
{
  // Check boundary condition options
  bc_options_.RequireBlockTypeIs(ParameterBlockType::ARRAY);

  // Get the initial density guess
  auto guess_params = params.GetParam("initial_guess");

  guess_params.RequireParameter("material_ids");
  material_ids_ = guess_params.GetParamVectorValue<int>("material_ids");

  guess_params.RequireParameter("values");
  densities_ = guess_params.GetParamVectorValue<double>("values");

  OpenSnLogicalErrorIf(material_ids_.size() != densities_.size(),
                       "The number of material IDs and values present in the initial guess "
                       "must be equivalent.");

  log.Log() << "\n***** inverse_solver Initialized *****\n";
}

void
InverseSolver::Initialize()
{
  solver_.Initialize();

  // Check density guess
  for (const auto& matid : material_ids_)
  {
    OpenSnLogicalErrorIf(solver_.GetMatID2XSMap().count(matid) == 0,
                         "Material ID " + std::to_string(matid) + " in the density map " +
                           "does not correspond to a material ID in the simulation.");
  }

  // Generate the measured signal to reconstruct
  SetForwardMode();
  Solve();
  ref_leakage_ = ComputeDetectorSignal();

  log.Log() << "\n********** Reference Signal **********\n"
               "Boundary         Value\n"
               "--------         -----\n";
  for (int i = 0; i < ref_leakage_.size(); ++i)
    log.Log() << std::left << std::setw(8) << detector_boundaries_[i] << "         " << std::left
              << std::setw(8) << std::setprecision(8) << ref_leakage_[i];
  log.Log() << "\n";
}

void
InverseSolver::Execute()
{
  // Set the densities to the initial guess.
  for (int i = 0; i < material_ids_.size(); ++i)
  {
    auto& xs = solver_.GetMatID2XSMap().at(material_ids_[i]);
    std::dynamic_pointer_cast<SingleStateMGXS>(xs)->SetScalingFactor(densities_[i]);
  }

  // Define the convergence normalization
  const auto norm = std::accumulate(ref_leakage_.begin(), ref_leakage_.end(), 0.0);

  // Start iterations
  double f, dphi, dphi_ell;
  double alpha, alpha_ell = alpha_max_;
  std::vector<double> df(densities_.size()), p(densities_.size());
  for (int nit = 0; nit < max_iterations_; ++nit)
  {
    // Compute objective function, its gradient, and step length
    f = EvaluateObjective();
    df = p = EvaluateGradient();

    if (line_search_)
    {
      std::transform(df.begin(), df.end(), p.begin(), std::negate{});
      dphi = std::inner_product(df.begin(), df.end(), p.begin(), 0.0);

      auto alpha0 = (nit == 0) ? alpha_max_ : alpha_ell * dphi_ell / dphi;
      alpha0 = std::min(alpha_max_, alpha0);
      alpha = BacktrackingLineSearch(alpha0, f, df, p);
    }
    else
    {
      alpha = alpha_max_;
    }
    
    // Update for next iteration
    densities_ = UpdateDensities(alpha, df);
    dphi_ell = dphi;
    alpha_ell = alpha;

    // Print iteration info
    if (nit == 0)
      log.Log() << "********** Iteration Summaries **********";

    std::stringstream ss;
    ss.precision(6);

    ss << "Iteration: " << std::setw(4) << std::left << nit << "    ";
    ss << "Obj Func: " << std::left << std::setw(12) << f << "    ";
    ss << "Conv:  " << std::left << std::setw(12) << f / norm << "    ";
    ss << "Alpha: " << std::left << std::setw(12) << alpha << "    ";
    ss << "Densities: [";
    for (int i = 0; i < densities_.size() - 1; ++i)
      ss << densities_[i] << " ";
    ss << densities_.back() << "]" << (f / norm < tolerance_ ? "  CONVERGED" : "");
    log.Log() << ss.str();

    if (f / norm < tolerance_)
      break;
  }
}

double
InverseSolver::EvaluateObjective() const
{
  // Compute the simulated leakage
  SetForwardMode();
  Solve();
  const auto leakage = ComputeDetectorSignal();

  // Evaluate the leakage mismatch
  double val = 0.0;
  for (int i = 0; i < ref_leakage_.size(); ++i)
    val += ref_leakage_[i] - leakage[i];
  return 0.5 * val * val;
}

std::vector<double>
InverseSolver::EvaluateGradient() const
{
  // Get the forward solutions
  const auto phi = solver_.PhiNewLocal();
  const auto psi = solver_.PsiNewLocal();

  // Solve the adjoint, get the adjoint solution
  SetAdjointMode();
  Solve();
  const auto psi_dagger = solver_.PsiNewLocal();

  // Reset to forward mode for later computations
  SetForwardMode();

  const auto& grid = solver_.Grid();
  const auto& discretization = solver_.SpatialDiscretization();
  const auto& unit_cell_matrices = solver_.GetUnitCellMatrices();
  const auto& transport_views = solver_.GetCellTransportViews();

  // Loop over groupsets
  int gs = 0;
  std::map<int, double> update;
  for (const auto& groupset : solver_.Groupsets())
  {
    const auto& uk_man = groupset.psi_uk_man_;

    const auto& quadrature = groupset.quadrature_;
    const auto& m2d = quadrature->GetMomentToDiscreteOperator();
    const auto& moment_map = quadrature->GetMomentToHarmonicsIndexMap();
    const auto& weights = quadrature->weights_;
    const auto num_gs_angles = quadrature->omegas_.size();

    const auto num_gs_groups = groupset.groups_.size();
    const auto first_gs_group = groupset.groups_.front().id_;

    const auto num_moments = solver_.NumMoments();

    // Loop over cells
    for (const auto& cell : grid.local_cells)
    {
      const auto& cell_mapping = discretization.GetCellMapping(cell);
      const auto& transport_view = transport_views[cell.local_id_];
      const auto& fe_values = unit_cell_matrices[cell.local_id_];
      const auto& xs = transport_view.XS();

      const auto matid = cell.material_id_;
      const auto density = xs.ScalingFactor();
      const auto& sigma_t = xs.SigmaTotal();
      const auto& S = xs.TransferMatrices();
      const auto& F = xs.ProductionMatrix();

      // Loop over cell nodes
      const auto num_cell_nodes = cell_mapping.NumNodes();
      for (int i = 0; i < num_cell_nodes; ++i)
      {
        const auto& V = fe_values.intV_shapeI[i];

        // Loop over groupset groups
        for (int gsg = 0; gsg < num_gs_groups; ++gsg)
        {
          const auto g = first_gs_group + gsg;
          const auto& sig_t = sigma_t[g];

          // Precompute source moments
          std::vector<double> q_moments(num_moments, 0.0);
          for (int m = 0; m < num_moments; ++m)
          {
            const auto ell = moment_map[m].ell;
            const auto uk_map = transport_view.MapDOF(i, m, 0);

            // Apply scattering
            double s = 0.0;
            if (ell < S.size())
              for (const auto& [_, gp, sig_ell] : S[ell].Row(g))
                s += sig_ell * phi[uk_map + gp];

            // Apply fission
            double f = 0.0;
            if (ell == 0 and xs.IsFissionable())
            {
              unsigned int gp = 0;
              for (const auto& sig_f : F[g])
                f += sig_f * phi[uk_map + gp];
            }

            q_moments[m] = s + f;

          } // for moment m

          // Apply the angular integration
          for (int n = 0; n < num_gs_angles; ++n)
          {
            const auto wt = weights[n] * V;
            const auto uk_map = discretization.MapDOF(cell, i, uk_man, n, gsg);

            // Compute discrete source from moment source
            double transfer = 0.0;
            for (unsigned int m = 0; m < solver_.NumMoments(); ++m)
              transfer += m2d[m][n] * q_moments[m];

            const double drho_A_psi = (sig_t * psi[gs][uk_map] - transfer) / density;
            update[matid] += wt * psi_dagger[gs][uk_map] * drho_A_psi;

          } // for angle n
        }   // for groupset group gsg
      }     // for node i
    }       // for cell

    ++gs;

  } // for groupset

  std::vector<double> local_outputs(material_ids_.size(), 0.0);
  for (int i = 0; i < material_ids_.size(); ++i)
    local_outputs[i] += update[material_ids_[i]];

  std::vector<double> global_outputs(material_ids_.size(), 0.0);
  mpi_comm.all_reduce(local_outputs.data(),
                      static_cast<int>(local_outputs.size()),
                      global_outputs.data(),
                      mpicpp_lite::op::sum<double>());
  return global_outputs;
}

double
InverseSolver::BacktrackingLineSearch(const double alpha0,
                                      const double f0,
                                      const std::vector<double>& df,
                                      const std::vector<double>& p) const
{
  // Define the stop criteria
  const auto m = std::inner_product(df.begin(), df.end(), p.begin(), 0.0);
  const auto t = -1.0e-4 * m;

  // Start the line search
  double alpha = alpha0;
  for (int j = 0; j < 5; ++j)
  {
    // Set the densities, then evaluate the objective functions
    UpdateDensities(alpha, df);
    const auto f = EvaluateObjective();

    // If sufficient decrease, exit the routine
    if (f0 - f >= alpha * t)
      break;

    // Otherwise, scale alpha down
    alpha *= 0.5;
  }
  return std::max(alpha, 0.001);
}

std::vector<double>
InverseSolver::ComputeDetectorSignal() const
{
  // Compute the leakage
  std::vector<uint64_t> bids;
  bids.reserve(detector_boundaries_.size());
  for (const auto& bname : detector_boundaries_)
    bids.emplace_back(boundary_map_.at(bname));
  const auto leakage_map = solver_.ComputeLeakage(bids);

  // Format the leakage
  std::vector<double> leakage;
  leakage.reserve(leakage_map.size());
  for (const auto& [bid, vals] : leakage_map)
    leakage.push_back(std::accumulate(vals.begin(), vals.end(), 0.0));
  return leakage;
}

std::vector<double>
InverseSolver::UpdateDensities(const double alpha, const std::vector<double>& drho) const
{
  std::vector<double> rho;
  for (int i = 0; i < material_ids_.size(); ++i)
  {
    // Apply step length to density
    auto update = alpha * drho[i];
    if (densities_[i] + update < 0.0)
      update = densities_[i] - 1.0e-6;
    if (std::fabs(update) > 0.1 * densities_[i])
      update = (update > 0.0 ? 0.1 : -0.1) * densities_[i];

    rho.push_back(densities_[i] + update);

    // Set material properties for new density
    auto& xs = solver_.GetMatID2XSMap().at(material_ids_[i]);
    std::dynamic_pointer_cast<SingleStateMGXS>(xs)->SetScalingFactor(rho[i]);
  }
  return rho;
}

void
InverseSolver::SetForwardMode() const
{
  ParameterBlock user_params;
  user_params.AddParameter("adjoint", false);
  user_params.AddParameter("clear_boundary_conditions", true);

  ParameterBlock bcs("boundary_conditions");
  for (int b = 0; b < bc_options_.NumParameters(); ++b)
    bcs.AddParameter(bc_options_.GetParam(b));
  bcs.ChangeToArray();
  user_params.AddParameter(bcs);

  auto params = LBSSolver::OptionsBlock();
  params.AssignParameters(user_params);
  solver_.SetOptions(params);
}

void
InverseSolver::SetAdjointMode() const
{
  const auto leakage = ComputeDetectorSignal();

  ParameterBlock user_params;
  user_params.AddParameter("adjoint", true);
  user_params.AddParameter("clear_boundary_conditions", true);

  ParameterBlock bcs("boundary_conditions");
  for (int i = 0; i < detector_boundaries_.size(); ++i)
  {
    ParameterBlock bc(std::to_string(i + 1));

    const auto bndry_name = detector_boundaries_[i];

    bc.AddParameter("name", bndry_name);
    bc.AddParameter("type", "isotropic");

    double mismatch = leakage[i] - ref_leakage_[i];
    std::vector<double> vals(solver_.NumGroups(), mismatch);
    ParameterBlock group_strength("group_strength", vals);
    bc.AddParameter(group_strength);

    bcs.AddParameter(bc);
  }
  bcs.ChangeToArray();
  user_params.AddParameter(bcs);

  auto params = LBSSolver::OptionsBlock();
  params.AssignParameters(user_params);
  solver_.SetOptions(params);
}

void
InverseSolver::Solve() const
{
  auto& ags_solver = *solver_.GetPrimaryAGSSolver();
  ags_solver.Setup();
  ags_solver.Solve();

  if (solver_.Options().use_precursors)
    solver_.ComputePrecursors();

  if (solver_.Options().adjoint)
    solver_.ReorientAdjointSolution();

  solver_.UpdateFieldFunctions();
}

} // namespace opensn::lbs