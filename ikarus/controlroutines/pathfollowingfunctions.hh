// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file pathfollowingfunctions.hh
 * \brief Defines structures and methods related to subsidiary functions for control routines.
 *
 * This file contains the declarations of the StandardArcLength, LoadControlWithSubsidiaryFunction,
 * and DisplacementControl structs, which represent subsidiary functions for arc-length, load control, and
 * displacement control methods, respectively. These functions are used in path-following control routines.
 */

#pragma once

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <Eigen/Core>

#include <ikarus/controlroutines/common.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/defaultfunctions.hh>

namespace Ikarus {

namespace Impl {

  template <typename NLS>
  void checkHasNoValidIDBCForceFunction() {
    static_assert(not Concepts::HasValidIDBCForceFunction<NLS>,
                  "The given path following function is not implemented to handle inhomogeneous Dirichlet BCs.");
  }

} // namespace Impl

/**
 * \struct SubsidiaryArgs
 * \brief Structure containing arguments for subsidiary functions.
 *
 * This structure holds various arguments used by subsidiary functions in control routines.
 * \ingroup  controlroutines
 */
struct SubsidiaryArgs
{
  double stepSize;              ///< The step size in the control routine.
  Eigen::VectorX<double> DD;    ///< The vector representing the solution increment.
  double Dlambda;               ///< The increment in the load factor.
  double f;                     ///< The value of the subsidiary function.
  Eigen::VectorX<double> dfdDD; ///< The derivative of the subsidiary function with respect to DD.
  double dfdDlambda;            ///< The derivative of the subsidiary function with respect to Dlambda.
  int currentStep;              ///< The current step index in the control routine.

  void setZero(const Concepts::EigenType auto& firstParameter) {
    stepSize = 0.0;
    DD.resizeLike(firstParameter);
    DD.setZero();
    Dlambda = 0.0;
    f       = 0.0;
    dfdDD.resizeLike(firstParameter);
    dfdDD.setZero();
    dfdDlambda  = 0.0;
    currentStep = 0;
  }
};

/**
 * \struct ArcLength
 * \brief Structure representing the subsidiary function for the standard arc-length method.
 *
 * This structure provides methods to evaluate the subsidiary function, perform initial prediction,
 * and perform intermediate prediction for the standard arc-length control method.
 *
 * \details The equation for the arc length method reads
 * \f[
 * f(\mathrm{D}\mathbf{D}, \mathrm{D} \lambda)=
 * \sqrt{||\mathrm{D}\mathbf{D}||^2+ \psi^2 (\mathrm{D} \lambda)^2 }- \hat{s}, \f]
 * where \f$\mathrm{D}\mathbf{D}\f$ is the increment of the solution vector and \f$\mathrm{D} \lambda\f$ is the load
 * factor increment. \f$\psi\f$ is the to-be-determined correction factor for the different dimensionalities between
 * \f$\mathrm{D}\mathbf{D}\f$ and \f$\mathrm{D} \lambda\f$. The scalar  \f$\hat{s} \f$ defines the requested size of
 * the step.
 * \ingroup  controlroutines
 */
struct ArcLength
{
  /**
   * \brief Evaluates the subsidiary function for the standard arc-length method.
   *
   * This method calculates the subsidiary function value and its derivatives for the given arguments and stores it in
   * the given args structure.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   */
  template <typename NLS>
  void operator()(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) const {
    const auto idbcDelta = idbcIncrement(req, nonlinearSolver, args.Dlambda);
    const auto root = sqrt(args.DD.squaredNorm() + idbcDelta.squaredNorm() + psi * psi * args.Dlambda * args.Dlambda);
    args.f          = root - args.stepSize;
    if (not computedInitialPredictor) {
      args.dfdDD.setZero();
      args.dfdDlambda = 0.0;
    } else {
      args.dfdDD       = args.DD / root;
      const auto Dhat0 = idbcDelta / args.Dlambda; // extracting Dhat0, assuming IDBC: Dhat = lambda * Dhat0
      args.dfdDlambda  = ((Dhat0.dot(Dhat0) + psi * psi) * args.Dlambda) / root;
    }
  }

  /**
   * \brief Performs the initial prediction for the standard arc-length method.
   *
   * \details This method initializes the prediction step for the standard arc-length method it computes \f$\psi\f$ and
   * computes initial \f$\mathrm{D}\mathbf{D}\f$ and \f$\mathrm{D} \lambda\f$.
   * Based on Eq. 4.13 of \cite rothAlgorithmenZurNichtlinearen2020d, the scaling factor psi is modified for better
   * convergence characteristics.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   * \ingroup  controlroutines
   */
  template <typename NLS>
  void initialPrediction(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) {
    auto& residual               = nonlinearSolver.residual();
    double dlambda               = 1.0;
    auto idbcDelta               = idbcIncrement(req, nonlinearSolver, dlambda);
    const auto Dhat0             = idbcDelta / dlambda;
    const auto idbcForceFunction = nonlinearSolver.idbcForceFunction();
    const auto updateFunction    = nonlinearSolver.updateFunction();
    req.parameter()              = dlambda;
    decltype(auto) R             = residual(req);
    decltype(auto) K             = derivative(residual)(req);
    auto linearSolver            = createSPDLinearSolverFromNonLinearSolver(nonlinearSolver);

    if constexpr (Concepts::HasValidIDBCForceFunction<NLS>)
      R += idbcForceFunction(req);

    linearSolver.analyzePattern(K);
    linearSolver.factorize(K);
    linearSolver.solve(args.DD, -R);

    const auto DD2    = args.DD.squaredNorm();
    const auto Dhat02 = Dhat0.squaredNorm();

    psi    = sqrt(DD2 + Dhat02);
    auto s = sqrt(DD2 + Dhat02 + psi * psi * dlambda * dlambda);

    args.DD      = args.DD * args.stepSize / s;
    args.Dlambda = args.stepSize / s;

    idbcDelta *= args.Dlambda / dlambda;

    updateFunction(req.globalSolution(), args.DD);
    req.parameter() = args.Dlambda;

    // modify the globalSolution() considering inhomogeneous Dirichlet BCs
    if constexpr (Concepts::HasValidIDBCForceFunction<NLS>)
      req.syncParameterAndGlobalSolution(updateFunction);

    // rescaling of psi
    psi = sqrt((idbcDelta.squaredNorm() + args.DD.squaredNorm()) / (args.Dlambda * args.Dlambda));

    computedInitialPredictor = true;
  }

  /**
   * \brief Performs intermediate prediction for the standard arc-length method.
   *
   * This method updates the prediction step for the standard arc-length method.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   */
  template <typename NLS>
  void intermediatePrediction(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) {
    if (not computedInitialPredictor)
      DUNE_THROW(Dune::InvalidStateException, "initialPrediction has to be called before intermediatePrediction.");
    nonlinearSolver.updateFunction()(req.globalSolution(), args.DD);
    req.parameter() += args.Dlambda;
    if constexpr (Concepts::HasValidIDBCForceFunction<NLS>)
      req.syncParameterAndGlobalSolution(nonlinearSolver.updateFunction());
  }

  /** \brief The name of the PathFollowing method. */
  constexpr std::string name() const { return "Arc length"; }

private:
  double psi{0.0};
  bool computedInitialPredictor{false};
};

/**
 * \struct LoadControlSubsidiaryFunction
 * \brief Structure representing the subsidiary function for the load control method.
 *
 * \details The equation for the load control method reads
 * \f[
 * f(\mathrm{D}\mathbf{D}, \mathrm{D} \lambda)=
 * \mathrm{D} \lambda -  \hat{s}, \f]
 * where \f$\mathrm{D}\mathbf{D}\f$ is the increment of the solution vector and \f$\mathrm{D} \lambda\f$ is the load
 * factor increment. The scalar   \f$\hat{s} \f$ defines the requested size of the step.
 * \ingroup  controlroutines
 */
struct LoadControlSubsidiaryFunction
{
  /**
   * \brief Evaluates the subsidiary function for the load control method.
   *
   * This method calculates the subsidiary function value and its derivatives for the given arguments.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   */
  template <typename NLS>
  void operator()(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) const {
    Impl::checkHasNoValidIDBCForceFunction<NLS>();
    args.f = args.Dlambda - args.stepSize;
    args.dfdDD.setZero();
    args.dfdDlambda = 1.0;
  }

  /**
   * \brief Performs initial prediction for the load control method.
   *
   * This method initializes the prediction step for the load control method.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   */
  template <typename NLS>
  void initialPrediction(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) {
    Impl::checkHasNoValidIDBCForceFunction<NLS>();
    args.Dlambda             = args.stepSize;
    req.parameter()          = args.Dlambda;
    computedInitialPredictor = true;
  }

  /**
   * \brief Performs intermediate prediction for the load control method.
   *
   * This method updates the prediction step for the load control method.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   */
  template <typename NLS>
  void intermediatePrediction(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) {
    Impl::checkHasNoValidIDBCForceFunction<NLS>();
    if (not computedInitialPredictor)
      DUNE_THROW(Dune::InvalidStateException, "initialPrediction has to be called before intermediatePrediction.");
    req.parameter() += args.Dlambda;
  }

  /** \brief The name of the PathFollowing method. */
  constexpr std::string name() const { return "Load Control"; }

private:
  bool computedInitialPredictor{false};
};

/**
 * \struct DisplacementControl
 * \brief Structure representing the subsidiary function for the displacement control method.
 *
 * \details The equation for the load control method reads
 * \f[
 * f(\mathrm{D}\mathbf{D}, \mathrm{D} \lambda)=
 * ||\mathrm{D}\mathbf{D}|| -  \hat{s}, \f]
 * where \f$\mathrm{D}\mathbf{D}\f$ is the increment of the solution vector and \f$\mathrm{D} \lambda\f$ is the load
 * factor increment. The scalar  \f$\hat{s} \f$ defines the requested size of the step.
 * \ingroup  controlroutines
 */
struct DisplacementControl
{
  /**
   * \brief Constructor for DisplacementControl.
   *
   * \param p_controlledIndices Vector containing the indices of the controlled degrees of freedom.
   */
  explicit DisplacementControl(std::vector<int> p_controlledIndices)
      : controlledIndices{std::move(p_controlledIndices)} {}

  /**
   * \brief Evaluates the subsidiary function for the displacement control method.
   *
   * This method calculates the subsidiary function value and its derivatives for the given arguments.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   */
  template <typename NLS>
  void operator()(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) const {
    Impl::checkHasNoValidIDBCForceFunction<NLS>();
    const auto controlledDOFsNorm = args.DD(controlledIndices).norm();
    args.f                        = controlledDOFsNorm - args.stepSize;
    args.dfdDlambda               = 0.0;
    args.dfdDD.setZero();
    args.dfdDD(controlledIndices) = args.DD(controlledIndices) / controlledDOFsNorm;
  }

  /**
   * \brief Performs initial prediction for the displacement control method.
   *
   * This method initializes the prediction step for the displacement control method.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   */
  template <typename NLS>
  void initialPrediction(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) {
    Impl::checkHasNoValidIDBCForceFunction<NLS>();
    args.DD(controlledIndices).array() = args.stepSize;
    req.globalSolution()               = args.DD;
    computedInitialPredictor           = true;
  }

  /**
   * \brief Performs intermediate prediction for the displacement control method.
   *
   * This method updates the prediction step for the displacement control method.
   *
   * \tparam NLS Type of the nonlinear solver.
   *
   * \param req The solution.
   * \param nonlinearSolver The nonlinear solver.
   * \param args The subsidiary function arguments.
   */
  template <typename NLS>
  void intermediatePrediction(typename NLS::Domain& req, NLS& nonlinearSolver, SubsidiaryArgs& args) {
    Impl::checkHasNoValidIDBCForceFunction<NLS>();
    if (not computedInitialPredictor)
      DUNE_THROW(Dune::InvalidStateException, "initialPrediction has to be called before intermediatePrediction.");
    req.globalSolution() += args.DD;
  }

  /** \brief The name of the PathFollowing method. */
  constexpr std::string name() const { return "Displacement Control"; }

private:
  std::vector<int> controlledIndices; /**< Vector containing the indices of the controlled degrees of freedom. */
  bool computedInitialPredictor{false};
};
} // namespace Ikarus
