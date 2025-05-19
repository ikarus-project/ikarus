// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverlogger.hh
 * \brief Listener implementation for logging non-linear solvers
 * \ingroup observer
 */

#pragma once
#include <spdlog/spdlog.h>

#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/listener/listener.hh>

namespace Ikarus {

/**
 * \brief Implementation of an observer for logging non-linear solvers.
 */
class NonLinearSolverLogger : public Listener
{
public:
  template <typename BC>
  NonLinearSolverLogger& subscribeTo(BC& bc) {
    this->subscribe(bc, [&](NonLinearSolverMessages message, const BC::State& state) { this->update(message, state); });
    return *this;
  }

  /**
   * \brief Implementation of the update method for logging control messages with a control routine state.
   *
   * \param message The received nonlinear solver message.
   * \param state The received nonlinear solver state.
   */
  void update(NonLinearSolverMessages message, const Concepts::NonLinearSolverState auto& state) {
    constexpr bool isDomainAVector = Concepts::EigenVector<typename std::remove_cvref_t<decltype(state)>::Domain>;
    switch (message) {
      case NonLinearSolverMessages::INIT:
        init<isDomainAVector>();
        break;
      case NonLinearSolverMessages::ITERATION_ENDED:
        iterationEnded<isDomainAVector>();
        break;
      case NonLinearSolverMessages::RESIDUALNORM_UPDATED:
        rNorm_ = state.information.residualNorm;
        break;
      case NonLinearSolverMessages::SOLUTION_CHANGED:
        if constexpr (not isDomainAVector)
          lambda_ = state.domain.parameter();
        break;
      case NonLinearSolverMessages::CORRECTIONNORM_UPDATED:
        dNorm_ = state.information.correctionNorm;
        break;
      case NonLinearSolverMessages::FINISHED_SUCESSFULLY:
        finishedSuccessfully(state.information.iterations);
        break;
      default:
        break;
    }
  }

private:
  int iters_{0};
  double dNorm_{0};
  double rNorm_{0};
  double lambda_{0};

  template <bool isDomainAVector>
  void init() {
    iters_ = 1;
    rNorm_ = 0.0;
    dNorm_ = 0.0;
    spdlog::info("Non-linear solver started:");
    if (not isDomainAVector) {
      spdlog::info("{:<11} {:<20} {:<20} {:<20}", "Ite", "normR", "normD", "lambda");
      spdlog::info("-------------------------------------------------------------------------------");
    } else {
      spdlog::info("{:<11} {:<20} {:<20}", "Ite", "normR", "normD");
      spdlog::info("-------------------------------------------------");
    }
  }

  template <bool isDomainAVector>
  void iterationEnded() {
    if (not isDomainAVector)
      spdlog::info("{} {:<10d} {:<20.2e} {:<20.2e} {:<20.2e}", "", iters_, rNorm_, dNorm_, lambda_);
    else
      spdlog::info("{} {:<10d} {:<20.2e} {:<20.2e}", "", iters_, rNorm_, dNorm_);
    ++iters_;
  }

  void finishedSuccessfully(int numberOfIterations) { spdlog::info("Number of iterations: {}", numberOfIterations); }
};
} // namespace Ikarus
