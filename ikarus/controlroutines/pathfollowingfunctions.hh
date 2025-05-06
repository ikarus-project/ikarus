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
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <Eigen/Core>

#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {
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
   * \param args The subsidiary function arguments.
   */
  void operator()(SubsidiaryArgs& args) const {
    if (psi) {
      const auto root = sqrt(args.DD.squaredNorm() + psi.value() * psi.value() * args.Dlambda * args.Dlambda);
      args.f          = root - args.stepSize;
      args.dfdDD      = args.DD / root;
      args.dfdDlambda = (psi.value() * psi.value() * args.Dlambda) / root;
    } else {
      DUNE_THROW(Dune::InvalidStateException, "You have to call initialPrediction first. Otherwise psi is not defined");
    }
  }

  /**
   * \brief Performs the initial prediction for the standard arc-length method.
   *
   * This method initializes the prediction step for the standard arc-length method it computes \f$\psi\f$ and
   * computes initial \f$\mathrm{D}\mathbf{D}\f$ and \f$\mathrm{D} \lambda\f$.
   *
   * \tparam F Type of the residual function.
   * \param residual The residual function.
   * \param args The subsidiary function arguments.
   * \param req The solution.
   * \ingroup  controlroutines
   */
  template <typename F>
  requires(Ikarus::Concepts::LinearSolverCheck<Ikarus::LinearSolver, typename F::Traits::template Range<1>,
                                               typename F::Domain::SolutionVectorType>)
  void initialPrediction(typename F::Domain& req, F& residual, SubsidiaryArgs& args) {
    SolverTypeTag solverTag;
    using JacobianType = std::remove_cvref_t<typename F::Traits::template Range<1>>;
    static_assert((traits::isSpecializationTypeAndNonTypes<Eigen::Matrix, JacobianType>::value) or
                      (traits::isSpecializationTypeNonTypeAndType<Eigen::SparseMatrix, JacobianType>::value),
                  "Linear solver not implemented for the chosen derivative type of the non-linear operator");

    if constexpr (traits::isSpecializationTypeAndNonTypes<Eigen::Matrix, JacobianType>::value)
      solverTag = SolverTypeTag::d_LDLT;
    else
      solverTag = SolverTypeTag::sd_SimplicialLDLT;

    req.parameter()  = 1.0;
    decltype(auto) R = residual(req);
    decltype(auto) K = derivative(residual)(req);

    auto linearSolver = LinearSolver(solverTag); // for the linear predictor step
    linearSolver.analyzePattern(K);
    linearSolver.factorize(K);
    linearSolver.solve(args.DD, -R);

    const auto DD2 = args.DD.squaredNorm();

    psi    = sqrt(DD2);
    auto s = sqrt(psi.value() * psi.value() + DD2);

    args.DD      = args.DD * args.stepSize / s;
    args.Dlambda = args.stepSize / s;

    req.globalSolution() = args.DD;
    req.parameter()      = args.Dlambda;
  }

  /**
   * \brief Performs intermediate prediction for the standard arc-length method.
   *
   * This method updates the prediction step for the standard arc-length method.
   *
   * \tparam F Type of the residual function.
   * \param residual The residual function.
   * \param args The subsidiary function arguments.
   * \param req The solution.
   */
  template <typename F>
  void intermediatePrediction(typename F::Domain& req, F& residual, SubsidiaryArgs& args) {
    req.globalSolution() += args.DD;
    req.parameter() += args.Dlambda;
  }

  /** \brief The name of the PathFollowing method. */
  constexpr auto name() const { return std::string("Arc length"); }

private:
  std::optional<double> psi;
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
   * \param args The subsidiary function arguments.
   */
  void operator()(SubsidiaryArgs& args) const {
    args.f = args.Dlambda - args.stepSize;
    args.dfdDD.setZero();
    args.dfdDlambda = 1.0;
  }

  /**
   * \brief Performs initial prediction for the load control method.
   *
   * This method initializes the prediction step for the load control method.
   *
   * \tparam F Type of the residual function.
   * \param residual The residual function.
   * \param args The subsidiary function arguments.
   * \param req The solution.
   */
  template <typename F>
  void initialPrediction(typename F::Domain& req, F& residual, SubsidiaryArgs& args) {
    args.Dlambda    = args.stepSize;
    req.parameter() = args.Dlambda;
  }

  /**
   * \brief Performs intermediate prediction for the load control method.
   *
   * This method updates the prediction step for the load control method.
   *
   * \tparam F Type of the residual function.
   * \param residual The residual function.
   * \param args The subsidiary function arguments.
   * \param req The solution.
   */
  template <typename F>
  void intermediatePrediction(typename F::Domain& req, F& residual, SubsidiaryArgs& args) {
    req.parameter() += args.Dlambda;
  }

  /** \brief The name of the PathFollowing method. */
  constexpr auto name() const { return std::string("Load Control"); }
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
   * \param args The subsidiary function arguments.
   */
  void operator()(SubsidiaryArgs& args) const {
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
   * \tparam F Type of the residual function.
   * \param residual The residual function.
   * \param args The subsidiary function arguments.
   * \param req The solution.
   */
  template <typename F>
  void initialPrediction(typename F::Domain& req, F& residual, SubsidiaryArgs& args) {
    args.DD(controlledIndices).array() = args.stepSize;
    req.globalSolution()               = args.DD;
  }

  /**
   * \brief Performs intermediate prediction for the displacement control method.
   *
   * This method updates the prediction step for the displacement control method.
   *
   * \tparam F Type of the residual function.
   * \param residual The residual function.
   * \param args The subsidiary function arguments.
   * \param req The solution.
   */
  template <typename F>
  void intermediatePrediction(typename F::Domain& req, F& residual, SubsidiaryArgs& args) {
    req.globalSolution() += args.DD;
  }

  /** \brief The name of the PathFollowing method. */
  constexpr auto name() const { return std::string("Displacement Control"); }

private:
  std::vector<int> controlledIndices; /**< Vector containing the indices of the controlled degrees of freedom. */
};
} // namespace Ikarus
