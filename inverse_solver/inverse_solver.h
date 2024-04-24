#pragma once

#include "opensn/modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include <vector>
#include <map>
#include <string>

namespace opensn::lbs
{

class InverseSolver : public opensn::Solver
{
public:
  static InputParameters GetInputParameters();
  explicit InverseSolver(const InputParameters& params);

  virtual ~InverseSolver() = default;

  /**
   * Initialize the inverse solver by initializing the underlying solver object,
   * performing checks, and computing the reference density signal.
   */
  void Initialize() override;
  void Execute() override;

  /**
   * Evaluate the objective function for a given computed leakage.
   */
  double EvaluateObjective() const;

  /**
   * Computes the density update by computing the density derivative of the
   * Lagrangian minimization functional.
   */
  std::vector<double> EvaluateGradient() const;

private:
  double BacktrackingLineSearch(double alpha0,
                                double f0,
                                const std::vector<double>& df,
                                const std::vector<double>& p) const;

  std::vector<double> ComputeDetectorSignal() const;

  /**
   * Sets material densities based on the currently set step length and the
   * density gradient vector.
   */
  std::vector<double> UpdateDensities(const double alpha, const std::vector<double>& drho) const;

  void SetForwardMode() const;
  void SetAdjointMode() const;

  /**
   * A generic steady state execution routine. This is equivalent to the
   * steady state executioner in LinearBoltzmannSolvers.
   */
  void Solve() const;

private:
  DiscreteOrdinatesSolver& solver_;
  const ParameterBlock bc_options_;

  const std::vector<std::string> detector_boundaries_;
  std::vector<double> ref_leakage_;
  std::vector<double> leakage_;

  std::vector<int> material_ids_;
  std::vector<double> densities_;

  const double alpha_max_;

  const unsigned int max_iterations_;
  const double tolerance_;

  const bool line_search_;
};

} // namespace opensn::lbs
