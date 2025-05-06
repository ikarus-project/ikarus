// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file loadcontrol.hh
 * \brief Defines the LoadControl class.
 */

#pragma once

#include <memory>

#include <dune/common/hybridutilities.hh>

#include <ikarus/controlroutines/controlinfos.hh>
#include <ikarus/controlroutines/controlroutinebase.hh>
#include <ikarus/controlroutines/controlroutinefactory.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>

namespace Ikarus {

template <typename NLS>
class LoadControl;

/**
 * \struct LoadControlConfig
 * \brief Config for the Load-Control control routine
 */
struct LoadControlConfig
{
  int loadSteps{};
  double tbegin{};
  double tEnd{};
};

/**
 * \brief Function to create a load control instance
 *
 * \tparam NLS Type of the nonlinear solver
 * \param config the provided config for the load control
 * \param nonlinearSolver the provided nonlinearsolver
 * \return LoadControl the newly created load control instance
 */
template <typename NLS>
auto createControlRoutine(const LoadControlConfig& config, NLS&& nonlinearSolver) {
  return LoadControl(std::forward<NLS>(nonlinearSolver), config.loadSteps, std::array{config.tbegin, config.tEnd});
}

/**
 * \class LoadControl
 * \brief The LoadControl control routine increases the parameter of the fe requirements given in run function and
 * solves the corresponding differentiable function f for its root and calls a nonlinear solver. \ingroup
 * controlroutines This class represents the LoadControl control routine. It increments the parameter of the fe
 * requirement and utilizes a nonlinear solver, such as Newton's method, to solve the resulting system at each step.
 *
 * \tparam NLS Type of the nonlinear solver used in the control routine.
 */
template <typename NLS>
class LoadControl : public ControlRoutineBase<typename NLS::DifferentiableFunction>
{
public:
  /** \brief The name of the LoadControl method. */
  constexpr auto name() const { return std::string("Load Control Method"); }

  /**
   * \brief Constructor for LoadControl.
   *
   * \param nonLinearSolver Shared pointer to the nonlinear solver.
   * \param loadSteps Number of load steps in the control routine.
   * \param tbeginEnd Array representing the range of load parameters [tbegin, tend].
   */
  LoadControl(const std::shared_ptr<NLS>& nonLinearSolver, int loadSteps, const std::array<double, 2>& tbeginEnd)
      : nonLinearSolver_{nonLinearSolver},
        loadSteps_{loadSteps},
        parameterBegin_{tbeginEnd[0]},
        parameterEnd_{tbeginEnd[1]},
        stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {}

  /**
   * \brief Executes the LoadControl routine.
   * \param x The solution.
   * \return ControlInformation structure containing information about the control results.
   */
  ControlInformation run(typename NLS::Domain& x);

  /* \brief returns the nonlinear solver */
  NLS& nonLinearSolver() { return *nonLinearSolver_; }

private:
  std::shared_ptr<NLS> nonLinearSolver_;
  int loadSteps_;
  double parameterBegin_;
  double parameterEnd_;
  double stepSize_;
};

} // namespace Ikarus

#include <ikarus/controlroutines/loadcontrol.inl>
