#pragma once

#include "opensn/framework/physics/solver_base/solver.h"
#include <petsctao.h>
#include <vector>
#include <map>
#include <string>

namespace opensn::lbs
{
class DiscreteOrdinatesSolver;

/**
 * An inverse problem solver for density reconstruction.
 */
class InverseSolver : public opensn::Solver
{
public:
  static InputParameters GetInputParameters();
  explicit InverseSolver(const InputParameters& params);

  ~InverseSolver() override;

  void Initialize() override;
  void Execute() override;

  PetscErrorCode EvaluateObjective(Vec x, PetscReal* f);
  PetscErrorCode EvaluateObjectiveAndGradient(Vec x, PetscReal* f, Vec df);

protected:
  void ExecuteSimple();
  void ExecuteTao();

  std::vector<double> ComputeDetectorSignal() const;
  void ComputeInnerProduct(const std::vector<double>& phi_fwd,
                           const std::vector<std::vector<double>>& psi_fwd,
                           const std::vector<std::vector<double>>& psi_adj,
                           Vec df) const;

  double BackTrackingLineSearch(Vec x, double alpha0, Vec df, double f0);
  void SetDensities(Vec x);

  InputParameters GetForwardOptions() const;
  InputParameters GetAdjointOptions(const std::vector<PetscReal>& leakage) const;
  void ExecuteSteadyState();

protected:
  Tao tao_;
  Vec solution_;
  Vec gradient_;

  std::vector<int> material_ids_;
  std::vector<uint64_t> detector_bids_;
  std::vector<double> ref_leakage_;

  DiscreteOrdinatesSolver& solver_;

  const ParameterBlock forward_bcs_;
  const std::vector<std::string> detector_bndrys_;

  const unsigned max_its_;
  const double tol_;
  const double alpha_;
  const bool line_search_;
  const unsigned max_ls_its_;
  const bool use_tao_;

  unsigned num_transport_solves_;

private:
  static std::string VecString(Vec x);
};

} // namespace opensn::lbs
