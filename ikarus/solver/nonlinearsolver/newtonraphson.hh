// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file newtonraphson.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverbase.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverstate.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

template <typename F, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
class NewtonRaphson;

struct NRSettings
{
  double tol{1e-8};
  int maxIter{20};
  int minIter{0};
};

/**
 * \struct NewtonRaphsonConfig
 * \brief Config for the Newton-Raphson solver.
 */
template <typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
struct NewtonRaphsonConfig
{
  using LinearSolver   = LS;
  using UpdateFunction = UF;
  NRSettings parameters{};
  LS linearSolver{};
  UF updateFunction{};

  template <typename UF2>
  auto rebindUpdateFunction(UF2&& updateFunction) const {
    NewtonRaphsonConfig<LS, UF2> settings{
        .parameters = parameters, .linearSolver = linearSolver, .updateFunction = std::forward<UF2>(updateFunction)};
    return settings;
  }

  template <typename F>
  using Solver = NewtonRaphson<F, LS, UF>;
};

// THE CTAD is broken for designated initializers in clang 16, when we drop support this can be simplified
#ifndef DOXYGEN
NewtonRaphsonConfig() -> NewtonRaphsonConfig<utils::SolverDefault, utils::UpdateDefault>;
NewtonRaphsonConfig(NRSettings) -> NewtonRaphsonConfig<utils::SolverDefault, utils::UpdateDefault>;

template <typename LS>
NewtonRaphsonConfig(NRSettings, LS) -> NewtonRaphsonConfig<LS, utils::UpdateDefault>;

template <typename LS, typename UF>
NewtonRaphsonConfig(NRSettings, LS, UF) -> NewtonRaphsonConfig<LS, UF>;

template <typename UF>
NewtonRaphsonConfig(NRSettings, utils::SolverDefault, UF) -> NewtonRaphsonConfig<utils::SolverDefault, UF>;
#endif

/**
 * \brief Function to create a NewtonRaphson solver instance.
 * \tparam F Type of the differentiable function to solve.
 * \tparam NRConfig Type of the nonlinear solver config.
 * \param config Config for the solver.
 * \param f Function to solve.
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename F, typename NRConfig>
requires traits::isSpecialization<NewtonRaphsonConfig, std::remove_cvref_t<NRConfig>>::value
auto createNonlinearSolver(NRConfig&& config, F&& f) {
  using LS           = std::remove_cvref_t<NRConfig>::LinearSolver;
  using UF           = std::remove_cvref_t<NRConfig>::UpdateFunction;
  auto solverFactory = []<class F2, class LS2, class UF2>(F2&& f2, LS2&& ls, UF2&& uf) {
    return std::make_shared<NewtonRaphson<std::remove_cvref_t<F2>, std::remove_cvref_t<LS2>, std::remove_cvref_t<UF2>>>(
        f2, std::forward<LS2>(ls), std::forward<UF2>(uf));
  };

  if constexpr (std::remove_cvref_t<F>::nDerivatives == 2) {
    auto solver = solverFactory(derivative(f), std::forward<NRConfig>(config).linearSolver,
                                std::forward<NRConfig>(config).updateFunction);
    solver->setup(config.parameters);
    return solver;
  } else {
    static_assert(std::remove_cvref_t<F>::nDerivatives >= 1,
                  "The number of derivatives in the differentiable function have to be more than 0");
    auto solver =
        solverFactory(f, std::forward<NRConfig>(config).linearSolver, std::forward<NRConfig>(config).updateFunction);
    ;

    solver->setup(std::forward<NRConfig>(config).parameters);
    return solver;
  }
}

/**
 * \class NewtonRaphson
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 * \tparam F Type of the differentiable function to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \relates makeNewtonRaphson
 * \ingroup solvers
 */
template <typename F, typename LS, typename UF>
class NewtonRaphson : public NonlinearSolverBase<F>
{
public:
  using Settings        = NRSettings;
  using SignatureTraits = typename F::Traits;
  using CorrectionType  = typename SignatureTraits::template Range<0>; ///< Type of the correction of x += deltaX.
  using JacobianType    = typename SignatureTraits::template Range<1>;
  ///< Compile-time boolean indicating if the linear solver satisfies the non-linear solver concept
  static constexpr bool isLinearSolver = Ikarus::Concepts::LinearSolverCheck<LS, JacobianType, CorrectionType>;

  using Domain = typename SignatureTraits::Domain; ///< Type representing the parameter vector of the function.

  using UpdateFunction         = UF; ///< Type representing the update function.
  using DifferentiableFunction = F;  ///< Type of the non-linear operator

  /**
   * \brief Constructor for NewtonRaphson.
   * \param residual residual to solve.
   * \param linearSolver Linear solver used internally (default is SolverDefault).
   * \param updateFunction Update function (default is UpdateDefault).
   */
  template <typename LS2 = LS, typename UF2 = UF>
  explicit NewtonRaphson(const DifferentiableFunction& residual, LS2&& linearSolver = {}, UF2&& updateFunction = {})
      : residualFunction_{residual},
        jacobianFunction_{derivative(residualFunction_)},
        linearSolver_{std::forward<LS2>(linearSolver)},
        updateFunction_{std::forward<UF2>(updateFunction)} {}

  /**
   * \brief Set up the solver with the given settings.
   * \param settings Newton-Raphson settings.
   */
  void setup(const Settings& settings) {
    if (settings.minIter > settings.maxIter)
      DUNE_THROW(Dune::InvalidStateException,
                 "Minimum number of iterations cannot be greater than maximum number of iterations");
    settings_ = settings;
  }

  /**
   * \brief Solve the nonlinear system.
   * \param x Where the solution should be stored.
   * \return Information about the solution process.
   */
  [[nodiscard(
      "The solve method returns information of the solution process. You should store this information and check if "
      "it was successful")]] Ikarus::NonLinearSolverInformation
  solve(Domain& x) {
    using enum NonLinearSolverMessages;
    this->notify(INIT);

    Ikarus::NonLinearSolverInformation solverInformation;
    solverInformation.success = true;

    decltype(auto) rx = residualFunction_(x);
    decltype(auto) Ax = jacobianFunction_(x);
    auto rNorm        = norm(rx);
    decltype(rNorm) dNorm;
    int iter{0};
    if constexpr (isLinearSolver)
      linearSolver_.analyzePattern(Ax);

    auto solverState = typename NewtonRaphson::State{.domain = x, .correction = correction_};

    while ((rNorm > settings_.tol && iter < settings_.maxIter) or iter < settings_.minIter) {
      this->notify(ITERATION_STARTED);
      if constexpr (isLinearSolver) {
        linearSolver_.factorize(Ax);
        linearSolver_.solve(correction_, -rx);
        dNorm = correction_.norm();
      } else {
        correction_ = -linearSolver_(rx, Ax);
        dNorm       = norm(correction_);
      }
      this->notify(CORRECTION_UPDATED, solverState);
      updateFunction_(x, correction_);
      this->notify(CORRECTIONNORM_UPDATED, static_cast<double>(dNorm));
      this->notify(SOLUTION_CHANGED);
      rx    = residualFunction_(x);
      Ax    = jacobianFunction_(x);
      rNorm = norm(rx);
      this->notify(RESIDUALNORM_UPDATED, static_cast<double>(rNorm));
      this->notify(ITERATION_ENDED);
      ++iter;
    }
    if (iter == settings_.maxIter)
      solverInformation.success = false;
    solverInformation.iterations     = iter;
    solverInformation.residualNorm   = static_cast<double>(rNorm);
    solverInformation.correctionNorm = static_cast<double>(dNorm);
    if (solverInformation.success)
      this->notify(FINISHED_SUCESSFULLY, iter);
    return solverInformation;
  }

  /**
   * \brief Access the function.
   * \return Reference to the function.
   */
  auto& residual() { return residualFunction_; }

private:
  DifferentiableFunction residualFunction_;
  typename DifferentiableFunction::Derivative jacobianFunction_;
  std::remove_cvref_t<CorrectionType> correction_;
  LS linearSolver_;
  UpdateFunction updateFunction_;
  Settings settings_;
};

/**
 * \brief Function to create a NewtonRaphson solver instance.
 * \tparam F Type of the function to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \param f Function to solve.
 * \param linearSolver Linear solver used internally (default is SolverDefault).
 * \param updateFunction Update function (default is UpdateDefault).
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename F, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
auto makeNewtonRaphson(const F& f, LS&& linearSolver = {}, UF&& updateFunction = {}) {
  return std::make_shared<NewtonRaphson<F, LS, UF>>(f, std::forward<LS>(linearSolver), std::move(updateFunction));
}

template <typename F, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
NewtonRaphson(const F& f, LS&& linearSolver = {},
              UF&& updateFunction = {}) -> NewtonRaphson<F, std::remove_cvref_t<LS>, std::remove_cvref_t<UF>>;

} // namespace Ikarus
