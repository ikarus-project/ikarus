// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file newtonraphson.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

#include <type_traits>
#include <utility>

#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverbase.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverstate.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus {

template <typename F, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault,
          typename IDBCF = utils::IDBCForceDefault>
class NewtonRaphsonWithSubsidiaryFunction;

struct NewtonRaphsonWithSubsidiaryFunctionSettings
{
  double tol{1e-8};
  int maxIter{20};
};

/**
 * \struct NewtonRaphsonWithSubsidiaryFunctionConfig
 * \brief Settings for the Newton-Raphson solver with subsidiary function.
 */
template <typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault,
          typename IDBCF = utils::IDBCForceDefault>
struct NewtonRaphsonWithSubsidiaryFunctionConfig
{
  using LinearSolver      = LS;
  using UpdateFunction    = UF;
  using IDBCForceFunction = IDBCF;
  NewtonRaphsonWithSubsidiaryFunctionSettings parameters{};
  LS linearSolver{};
  UF updateFunction{};
  IDBCF idbcForceFunction{};

  template <typename UF2>
  auto rebindUpdateFunction(UF2&& updateFunction) const {
    NewtonRaphsonWithSubsidiaryFunctionConfig<LS, UF2, IDBCF> settings{
        .parameters        = parameters,
        .linearSolver      = linearSolver,
        .updateFunction    = std::forward<UF2>(updateFunction),
        .idbcForceFunction = idbcForceFunction};
    return settings;
  }

  template <typename IDBCF2>
  auto rebindIDBCForceFunction(IDBCF2&& idbcForceFunction) const {
    NewtonRaphsonWithSubsidiaryFunctionConfig<LS, UF, IDBCF2> settings{
        .parameters        = parameters,
        .linearSolver      = linearSolver,
        .updateFunction    = updateFunction,
        .idbcForceFunction = std::forward<IDBCF2>(idbcForceFunction)};
    return settings;
  }

  template <typename F>
  using Solver = NewtonRaphsonWithSubsidiaryFunction<F, LS, UF>;
};

// THE CTAD is broken for designated initializers in clang 16, when we drop support this can be simplified
#ifndef DOXYGEN
NewtonRaphsonWithSubsidiaryFunctionConfig()
    -> NewtonRaphsonWithSubsidiaryFunctionConfig<utils::SolverDefault, utils::UpdateDefault, utils::IDBCForceDefault>;

NewtonRaphsonWithSubsidiaryFunctionConfig(NewtonRaphsonWithSubsidiaryFunctionSettings)
    -> NewtonRaphsonWithSubsidiaryFunctionConfig<utils::SolverDefault, utils::UpdateDefault, utils::IDBCForceDefault>;

template <typename LS>
NewtonRaphsonWithSubsidiaryFunctionConfig(NewtonRaphsonWithSubsidiaryFunctionSettings, LS)
    -> NewtonRaphsonWithSubsidiaryFunctionConfig<LS, utils::UpdateDefault, utils::IDBCForceDefault>;

template <typename UF>
NewtonRaphsonWithSubsidiaryFunctionConfig(NewtonRaphsonWithSubsidiaryFunctionSettings, utils::SolverDefault, UF)
    -> NewtonRaphsonWithSubsidiaryFunctionConfig<utils::SolverDefault, UF, utils::IDBCForceDefault>;

template <typename IDBCF>
NewtonRaphsonWithSubsidiaryFunctionConfig(NewtonRaphsonWithSubsidiaryFunctionSettings, utils::SolverDefault,
                                          utils::UpdateDefault, IDBCF)
    -> NewtonRaphsonWithSubsidiaryFunctionConfig<utils::SolverDefault, utils::UpdateDefault, IDBCF>;

template <typename LS, typename UF>
NewtonRaphsonWithSubsidiaryFunctionConfig(NewtonRaphsonWithSubsidiaryFunctionSettings, LS, UF)
    -> NewtonRaphsonWithSubsidiaryFunctionConfig<LS, UF, utils::IDBCForceDefault>;

template <typename LS, typename UF, typename IDBCF>
NewtonRaphsonWithSubsidiaryFunctionConfig(NewtonRaphsonWithSubsidiaryFunctionSettings, LS, UF,
                                          IDBCF) -> NewtonRaphsonWithSubsidiaryFunctionConfig<LS, UF, IDBCF>;
#endif
/**
 * \brief Function to create a NewtonRaphson with subsidiary function solver instance.
 * \tparam F Type of the differentiable function to solve.
 * \tparam NRConfig Type of the nonlinear solver config.
 * \param config Config for the solver.
 * \param f Function to solve.
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename F, typename NRConfig>
requires traits::isSpecialization<NewtonRaphsonWithSubsidiaryFunctionConfig, std::remove_cvref_t<NRConfig>>::value
auto createNonlinearSolver(NRConfig&& config, F&& f) {
  using LS           = std::remove_cvref_t<NRConfig>::LinearSolver;
  using UF           = std::remove_cvref_t<NRConfig>::UpdateFunction;
  using IDBCF        = std::remove_cvref_t<NRConfig>::IDBCForceFunction;
  auto solverFactory = []<class F2, class LS2, class UF2, class IDBCF2>(F2&& f2, LS2&& ls, UF2&& uf, IDBCF2&& if_) {
    return std::make_shared<NewtonRaphsonWithSubsidiaryFunction<std::remove_cvref_t<F2>, std::remove_cvref_t<LS2>,
                                                                std::remove_cvref_t<UF2>, std::remove_cvref_t<IDBCF2>>>(
        f2, std::forward<LS2>(ls), std::forward<UF2>(uf), std::forward<IDBCF2>(if_));
  };

  if constexpr (std::remove_cvref_t<F>::nDerivatives == 2) {
    auto solver =
        solverFactory(derivative(f), std::forward<NRConfig>(config).linearSolver,
                      std::forward<NRConfig>(config).updateFunction, std::forward<NRConfig>(config).idbcForceFunction);
    solver->setup(config.parameters);
    return solver;
  } else {
    static_assert(std::remove_cvref_t<F>::nDerivatives >= 1,
                  "The number of derivatives in the differentiable function have to be more than 0");
    auto solver =
        solverFactory(f, std::forward<NRConfig>(config).linearSolver, std::forward<NRConfig>(config).updateFunction,
                      std::forward<NRConfig>(config).idbcForceFunction);

    solver->setup(std::forward<NRConfig>(config).parameters);
    return solver;
  }
}

/**
 * \brief Newton-Raphson solver with subsidiary function.
 *
 * This class provides a Newton-Raphson solver for solving nonlinear systems with a subsidiary function.
 * It uses a linear solver to handle the linear system arising in each iteration.
 *
 * \tparam F Type of the function.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \tparam IDBCF Type of the force function to handle inhomogeneous Dirichlet BCs (defaults to IDBCForceDefault).
 */
template <typename F, typename LS, typename UF, typename IDBCF>
class NewtonRaphsonWithSubsidiaryFunction : public NonlinearSolverBase<F>
{
public:
  using Settings        = NewtonRaphsonWithSubsidiaryFunctionSettings;
  using SignatureTraits = typename F::Traits;

  using Domain         = typename SignatureTraits::Domain; ///< Type representing the parameter vector of the Function.
  using CorrectionType = typename SignatureTraits::template Range<0>; ///< Type of the correction of x += deltaX.
  using JacobianType   = typename SignatureTraits::template Range<1>;
  ///< Compile-time boolean indicating if the linear solver satisfies the non-linear solver concept
  static constexpr bool isLinearSolver = Ikarus::Concepts::LinearSolverCheck<LS, JacobianType, CorrectionType>;

  ///< Type representing the update function.
  using UpdateFunction         = UF;
  using DifferentiableFunction = F;     ///< Type of the non-linear operator
  using IDBCForceFunction      = IDBCF; ///< Type representing the force function to handle inhomogeneous Dirichlet BCs.

  /**
   * \brief Constructor for NewtonRaphsonWithSubsidiaryFunction.
   * \param residual residual to solve.
   * \param linearSolver Linear solver used internally (default is SolverDefault).
   * \param updateFunction Update function (default is UpdateDefault).
   * \param idbcForceFunction Force function to handle inhomogeneous Dirichlet BCs (default is IDBCForceDefault).
   */
  template <typename LS2 = LS, typename UF2 = UF, typename IDBCF2 = IDBCF>
  explicit NewtonRaphsonWithSubsidiaryFunction(const DifferentiableFunction& residual, LS2&& linearSolver = {},
                                               UF2&& updateFunction = {}, IDBCF2&& idbcForceFunction = {})
      : residualFunction_{residual},
        jacobianFunction_{derivative(residualFunction_)},
        linearSolver_{std::forward<LS2>(linearSolver)},
        updateFunction_{std::forward<UF2>(updateFunction)},
        idbcForceFunction_{std::forward<IDBCF2>(idbcForceFunction)} {}

  /**
   * \brief Setup the Newton-Raphson solver with subsidiary function.
   *
   * \param p_settings Settings for the solver.
   */
  void setup(const Settings& settings) { settings_ = settings; }

  /**
   * \brief Solve the nonlinear system using the Newton-Raphson method with subsidiary function.
   *
   * \tparam SubsidiaryType Type of the subsidiary function.
   * \param subsidiaryFunction Subsidiary function to be solved.
   * \param req Where the solution should be stored.
   * \param subsidiaryArgs Additional arguments for the subsidiary function.
   * \param dxPredictor Predictor for the solution increment (default is NoPredictor).
   * \return Information about the solution process.
   */
  template <typename SubsidiaryType>
  [[nodiscard(
      "The solve method returns information of the solution process. You should store this information and check if "
      "it was successful")]] NonLinearSolverInformation
  solve(Domain& req, SubsidiaryType&& subsidiaryFunction, SubsidiaryArgs& subsidiaryArgs) {
    using enum NonLinearSolverMessages;

    Ikarus::NonLinearSolverInformation solverInformation{};
    auto state = typename NewtonRaphsonWithSubsidiaryFunction::State{
        .domain = req, .correction = correction_, .information = solverInformation};
    this->notify(INIT, state);

    solverInformation.success = true;

    auto& lambda = req.parameter();
    auto& x      = req.globalSolution();

    /// Determine Fext0
    /// It is assumed that Fext = Fext0 * lambda such that dRdlambda = Fext0
    /// Generalization for Fext0 = Fext0(lambda) is not implemented
    /// Forces due to inhomogeneous Dirichlet BCs is treated separately

    auto lambdaDummy = lambda;
    auto DDummy      = x;
    x.setZero();
    lambda            = 1.0;
    decltype(auto) rx = residualFunction_(req);
    const auto Fext0  = (-rx).eval(); // dRdlambda(lambda)
    lambda            = lambdaDummy;
    x                 = DDummy;

    rx                = residualFunction_(req);
    decltype(auto) Ax = jacobianFunction_(req);

    Eigen::VectorXd deltaD;
    deltaD.resizeLike(rx);
    deltaD.setZero();

    subsidiaryArgs.dfdDD.resizeLike(Fext0);

    subsidiaryFunction(req, *this, subsidiaryArgs);
    auto rNorm                     = sqrt(rx.dot(rx));
    solverInformation.residualNorm = static_cast<double>(rNorm);
    decltype(rNorm) dNorm;
    int iter{0};
    if constexpr (isLinearSolver)
      linearSolver_.analyzePattern(Ax);

    Eigen::MatrixX2<double> residual2d, sol2d;

    /// Iterative solving scheme
    while (rNorm > settings_.tol && iter < settings_.maxIter) {
      this->notify(ITERATION_STARTED, state);

      /// Two-step solving procedure
      residual2d.resize(rx.rows(), 2);
      sol2d.resize(rx.rows(), 2);
      if constexpr (not std::same_as<IDBCForceFunction, utils::IDBCForceDefault>)
        residual2d << -rx, Fext0 - (idbcForceFunction_(req));
      else
        residual2d << -rx, Fext0;

      if constexpr (isLinearSolver) {
        linearSolver_.factorize(Ax);
        linearSolver_.solve(sol2d, residual2d);
      } else {
        sol2d = linearSolver_(residual2d, Ax);
      }

      const auto& deltaDR = sol2d.col(0).eval();
      const auto& deltaDL = sol2d.col(1).eval();

      const double deltalambda = (-subsidiaryArgs.f - subsidiaryArgs.dfdDD.dot(deltaDR)) /
                                 (subsidiaryArgs.dfdDD.dot(deltaDL) + subsidiaryArgs.dfdDlambda);
      deltaD      = deltaDR + deltalambda * deltaDL;
      correction_ = deltaD;

      this->notify(CORRECTION_UPDATED, state);

      updateFunction_(x, deltaD);
      subsidiaryArgs.DD += deltaD;

      lambda += deltalambda;
      subsidiaryArgs.Dlambda += deltalambda;

      if constexpr (not std::same_as<IDBCForceFunction, utils::IDBCForceDefault>)
        req.syncParameterAndGlobalSolution(updateFunction_);

      dNorm                            = sqrt(deltaD.dot(deltaD) + deltalambda * deltalambda);
      solverInformation.correctionNorm = static_cast<double>(dNorm);

      rx = residualFunction_(req);
      Ax = jacobianFunction_(req);
      subsidiaryFunction(req, *this, subsidiaryArgs);

      rNorm                          = sqrt(rx.dot(rx) + subsidiaryArgs.f * subsidiaryArgs.f);
      solverInformation.residualNorm = static_cast<double>(rNorm);

      this->notify(SOLUTION_CHANGED, state);
      this->notify(CORRECTIONNORM_UPDATED, state);
      this->notify(RESIDUALNORM_UPDATED, state);

      ++iter;
      solverInformation.iterations = iter;
      this->notify(ITERATION_ENDED, state);
    }

    if (iter == settings_.maxIter)
      solverInformation.success = false;
    solverInformation.iterations = iter;
    if (solverInformation.success)
      this->notify(FINISHED_SUCESSFULLY, state);
    return solverInformation;
  }

  /**
   * \brief Access the residual function.
   * \return Reference to the residual function.
   */
  auto& residual() { return residualFunction_; }

  /**
   * \brief Access the function.
   * \return Reference to the function.
   */
  const auto& residual() const { return residualFunction_; }

  /**
   * \brief Access the update function.
   * \return Reference to the function.
   */
  const UpdateFunction& updateFunction() const { return updateFunction_; }

  /**
   * \brief Access the force function calculating internal forces due to inhomogeneous Dirichlet BCs.
   * \return Reference to the function.
   */
  const IDBCForceFunction& idbcForceFunction() const { return idbcForceFunction_; }

private:
  DifferentiableFunction residualFunction_;
  typename DifferentiableFunction::Derivative jacobianFunction_;
  std::remove_cvref_t<CorrectionType> correction_;
  LS linearSolver_;
  UpdateFunction updateFunction_;
  IDBCForceFunction idbcForceFunction_;
  Settings settings_;
};
/**
 * \brief Function to create a NewtonRaphson with subsidiary function solver instance.
 * \tparam F Type of the function to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \param f Function to solve.
 * \param linearSolver Linear solver used internally (default is SolverDefault).
 * \param updateFunction Update function (default is UpdateDefault).
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename F, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault,
          typename IDBCF = utils::IDBCForceDefault>
auto makeNewtonRaphsonWithSubsidiaryFunction(const F& f, LS&& linearSolver = {}, UF&& updateFunction = {},
                                             IDBCF&& idbcForceFunction = {}) {
  return std::make_shared<NewtonRaphsonWithSubsidiaryFunction<F, LS, UF, IDBCF>>(
      f, std::forward<LS>(linearSolver), std::move(updateFunction), std::move(idbcForceFunction));
}

template <typename F, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault,
          typename IDBCF = utils::IDBCForceDefault>
NewtonRaphsonWithSubsidiaryFunction(const F& f, LS&& linearSolver = {}, UF&& updateFunction = {},
                                    IDBCF&& idbcForceFunction = {})
    -> NewtonRaphsonWithSubsidiaryFunction<F, std::remove_cvref_t<LS>, std::remove_cvref_t<UF>,
                                           std::remove_cvref_t<IDBCF>>;

} // namespace Ikarus
