// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file loadcontrol.hh
 * \brief Defines the LoadControl class.
 */

#pragma once

#include <memory>

#include <ikarus/controlroutines/controlstate.hh>
#include <ikarus/solver/nonlinearsolver/solverstate.hh>
#include <ikarus/utils/observer/observable.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

/**
 * \class LoadControl
 * \brief The LoadControl control routine increases the last parameter of a nonlinear operator and calls a nonlinear
 * solver.
 *   \ingroup controlroutines
 * This class represents the LoadControl control routine. It increments the last parameter of a nonlinear operator
 * and utilizes a nonlinear solver, such as Newton's method, to solve the resulting system at each step.
 *
 * \tparam NLS Type of the nonlinear solver used in the control routine.
 */
template <typename NLS>
requires(Concepts::NonLinearSolverCheckForPathFollowing<NLS>)
class LoadControl : public ControlObservable
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
  LoadControl(const std::shared_ptr<NLS>& nonLinearSolver, size_t loadSteps, const std::array<double, 2>& tbeginEnd)
      : nonLinearSolver_{nonLinearSolver},
        loadSteps_{loadSteps},
        parameterBegin_{tbeginEnd[0]},
        parameterEnd_{tbeginEnd[1]},
        stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {
    static_assert(
        requires {
          nonLinearSolver_->nonLinearOperator().lastParameter() = 0.0;
          nonLinearSolver_->nonLinearOperator().lastParameter() += 0.0;
        }, "The last parameter (load factor) must be assignable and incrementable with a double!");
    if (loadSteps_ == 0)
      DUNE_THROW(Dune::InvalidStateException, "Number of load steps should be greater than zero.");
  }

  /**
   * \brief Executes the LoadControl routine.
   *
   * \return ControlState structure containing information about the control results.
   */
  ControlState run();

  /* \brief returns the nonlinear solver */
  NLS& nonlinearSolver() { return *nonLinearSolver_; }

private:
  std::shared_ptr<NLS> nonLinearSolver_;
  size_t loadSteps_;
  double parameterBegin_;
  double parameterEnd_;
  double stepSize_;

  /**
   * \brief A wrapper function to update controlState after a nonlinear solver is executed successfully.
   * The resulting controlState is then passed for further notifications.
   */
  void updateAndNotifyControlState(ControlState& controlState, typename NLS::NonLinearOperator& nonOp,
                                   const NonLinearSolverState& solverState) {
    controlState.solverStates.push_back(solverState);
    controlState.sol    = &nonOp.firstParameter();
    controlState.lambda = nonOp.lastParameter();
    this->notify(ControlMessages::SOLUTION_CHANGED, controlState);
    this->notify(ControlMessages::STEP_ENDED, controlState);
  }
};

} // namespace Ikarus

#include <ikarus/controlroutines/loadcontrol.inl>
