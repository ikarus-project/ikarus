// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file trustregion.hh
 * \brief Implementation of the Trust-region method for solving nonlinear equations.
 */

#pragma once

#include <iosfwd>
#include <type_traits>

#include <dune/common/float_cmp.hh>
#include <dune/functions/common/signature.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Sparse>

#include <ikarus/linearalgebra/truncatedconjugategradient.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverbase.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

/**
 * \enum PreConditioner
 * \brief Enumeration of available preconditioners for the trust region solver.
 */
enum class PreConditioner
{
  IncompleteCholesky,
  IdentityPreconditioner,
  DiagonalPreconditioner
};

struct TRSettings
{
  int verbosity    = 5;                                       ///< Verbosity level.
  double maxtime   = std::numeric_limits<double>::infinity(); ///< Maximum allowable time for solving.
  int minIter      = 3;                                       ///< Minimum number of iterations.
  int maxIter      = 1000;                                    ///< Maximum number of iterations.
  int debug        = 0;                                       ///< Debugging flag.
  double grad_tol  = 1e-6;                                    ///< Gradient tolerance.
  double corr_tol  = 1e-6;                                    ///< Correction tolerance.
  double rho_prime = 0.01;                                    ///< Rho prime value.
  bool useRand     = false;                                   ///< Flag for using random correction predictor.
  double rho_reg   = 1e6;                                     ///< Regularization value for rho.
  double Delta_bar = std::numeric_limits<double>::infinity(); ///< Maximum trust region radius.
  double Delta0    = 10;                                      ///< Initial trust region radius.
};

/**
 * \struct TrustRegionSettings
 * \brief Configuration settings for the TrustRegion solver.
 */

template <PreConditioner preConditioner = PreConditioner::IncompleteCholesky, typename UF = utils::UpdateDefault,
          typename IDBCF = utils::IDBCForceDefault>
struct TrustRegionConfig
{
  static_assert(std::copy_constructible<UF>,
                " The update function should be copy constructable. If it is a lambda wrap it in a std::function");

  using UpdateFunction    = UF;
  using IDBCForceFunction = IDBCF;

  TRSettings parameters;
  static constexpr PreConditioner preConditionerType = preConditioner;
  UF updateFunction;
  IDBCF idbcForceFunction{};
  template <typename UF2>
  auto rebindUpdateFunction(UF2&& updateFunction) const {
    TrustRegionConfig<preConditioner, UF2, IDBCF> settings{.parameters        = parameters,
                                                           .updateFunction    = std::forward<UF2>(updateFunction),
                                                           .idbcForceFunction = idbcForceFunction};
    return settings;
  }
  template <typename IDBC2>
  auto rebindIDBCForceFunction(IDBC2&& idbcForceFunction) const {
    TrustRegionConfig<preConditioner, UF, IDBC2> settings{.parameters        = parameters,
                                                          .updateFunction    = updateFunction,
                                                          .idbcForceFunction = std::forward<IDBC2>(idbcForceFunction)};
    return settings;
  }
};

template <typename F, PreConditioner preConditioner = PreConditioner::IncompleteCholesky,
          typename UF = utils::UpdateDefault, typename IDBCF = utils::IDBCForceDefault>
class TrustRegion;

/**
 * \brief Function to create a trust region non-linear solver
 * \tparam F Type of the differentiable function to solve.
 * \tparam TRConfig Type of the nonlinear solver config.
 * \param config Config for the solver.
 * \param f Function to solve.
 * \return Shared pointer to the TrustRegion solver instance.
 */
template <typename F, typename TRConfig>
requires traits::isSpecializationNonTypeAndTypes<TrustRegionConfig, std::remove_cvref_t<TRConfig>>::value
auto createNonlinearSolver(TRConfig&& config, F&& f) {
  static constexpr PreConditioner preConditioner = std::remove_cvref_t<TRConfig>::preConditionerType;
  using UF                                       = std::remove_cvref_t<TRConfig>::UpdateFunction;
  using IDBCF                                    = std::remove_cvref_t<TRConfig>::IDBCForceFunction;
  static_assert(std::remove_cvref_t<F>::nDerivatives == 2,
                "The number of derivatives in the DifferentiableFunction have to be exactly 2.");
  auto solver = std::make_shared<TrustRegion<F, preConditioner, UF, IDBCF>>(
      f, std::forward<TRConfig>(config).updateFunction, std::forward<TRConfig>(config).idbcForceFunction);

  solver->setup(config.parameters);
  return solver;
}

/**
 * \enum StopReason
 * \brief Enumeration of reasons for stopping the TrustRegion solver.
 */
enum class StopReason
{
  gradientNormTolReached,
  correctionNormTolReached,
  maximumTimeReached,
  maximumIterationsReached,
  dontStop
};

/**
 * \struct AlgoInfo
 * \brief Additional information about the TrustRegion algorithm.
 */
struct AlgoInfo
{
  int consecutive_TRplus   = 0;                 ///< Consecutive trust region increases.
  int consecutive_TRminus  = 0;                 ///< Consecutive trust region decreases.
  int consecutive_Rejected = 0;                 ///< Consecutive rejected proposals.
  std::string stopReasonString;                 ///< String describing the stopping reason.
  std::string trstr;                            ///< Trust region change string (TR+, TR-).
  std::string accstr;                           ///< Acceptance string (REJ, acc).
  std::string randomPredictionString;           ///< Random prediction string.
  std::string cauchystr = "                  "; ///< Used Cauchy step string.
  bool acceptProposal;                          ///< Flag indicating whether the proposal is accepted.
  bool used_cauchy = false;                     ///< Flag indicating whether Cauchy point was used.
  StopReason stop{StopReason::dontStop};        ///< Stopping reason.
  std::string reasonString;                     ///< String describing the stopping reason.
};

/**
 * \struct Stats
 * \brief Information about the TrustRegion solver.
 */
struct Stats
{
  double gradNorm{1};
  double etaNorm{1};
  double time{0};
  double energy{0};
  double energyProposal{0};
  double rho{0};
  int outerIter{0};
  int innerIterSum{0};
};

/**
* \class TrustRegion
* \brief Trust Region solver for non-linear optimization problems.
* \details Refer to \cite trustregion for details of the algorithm.

* This code is heavily inspired by the trust-region implementation of
<a href="https://github.com/NicolasBoumal/manopt/blob/master/manopt/solvers/trustregions/trustregions.m">Manopt</a>.
* \ingroup solvers
* \tparam F Type of the differentiable function to solve.
* \tparam preConditioner Type of preconditioner to use (default is IncompleteCholesky).
* \tparam UF Type of the update function
* \tparam IDBCF Type of the force function to handle inhomogeneous Dirichlet BCs (defaults to IDBCForceDefault).
*/
template <typename F, PreConditioner preConditioner, typename UF, typename IDBCF>
class TrustRegion : public NonlinearSolverBase<F>
{
public:
  using Settings = TRSettings; ///< Type of the settings for the TrustRegion solver

  using FTraits                = typename F::Traits;
  using Domain                 = typename FTraits::Domain;            ///< Type of the parameter vector of
  using CorrectionType         = typename FTraits::template Range<1>; ///< Type of the correction of x += deltaX.
  using UpdateFunction         = UF;                                  ///< Type of the update function.
  using DifferentiableFunction = F;                                   ///< Type of function to minimize
  using IDBCForceFunction      = IDBCF; ///< Type representing the force function to handle inhomogeneous Dirichlet BCs.

  using EnergyType   = typename FTraits::template Range<0>; ///< Type of the scalar cost
  using GradientType = typename FTraits::template Range<1>; ///< Type of the gradient vector
  using HessianType  = typename FTraits::template Range<2>; ///< Type of the Hessian matrix
  using JacobianType = HessianType;
  /**
   * \brief Constructs a TrustRegion solver instance.
   * \param f Function to solve.
   * \param updateFunction Update function
   * \param idbcForceFunction Force function to handle inhomogeneous Dirichlet BCs (default is IDBCForceDefault).
   */
  template <typename UF2 = UF, typename IDBCF2 = IDBCF>
  explicit TrustRegion(const F& f, UF2&& updateFunction = {}, IDBCF2&& idbcForceFunction = {})
      : energyFunction_{f},
        updateFunction_{std::forward<UF2>(updateFunction)},
        idbcForceFunction_{std::forward<IDBCF2>(idbcForceFunction)} {
    static_assert(std::same_as<IDBCForceFunction, utils::IDBCForceDefault>,
                  "Trust Region Method is not implemented to handle inhomogeneous Dirichlet BCs.");
  }

  /**
   * \brief Sets up the TrustRegion solver with the provided settings and checks feasibility.
   * \param p_settings TrustRegionSettings containing the solver configuration.
   */
  void setup(const Settings& settings) {
    settings_ = settings;
    assert(settings_.rho_prime < 0.25 && "options.rho_prime must be strictly smaller than 1/4.");
    assert(settings_.Delta_bar > 0 && "options.Delta_bar must be positive.");
    assert(settings_.Delta0 > 0 && settings_.Delta0 < settings_.Delta_bar &&
           "options.Delta0 must be positive and smaller than Delta_bar.");
  }

  /**
   * \brief Solves the nonlinear optimization problem using the TrustRegion algorithm.
   * \param x the solution.
   * \param stepSize the step size of the control routine (defaults to 0.0)
   * \return NonLinearSolverInformation containing information about the solver result.
   */
  [[nodiscard]] NonLinearSolverInformation solve(Domain& x, double stepSize = 0.0) {
    using enum NonLinearSolverMessages;

    NonLinearSolverInformation solverInformation;
    auto state = typename TrustRegion::State{.domain = x, .correction = eta_, .information = solverInformation};

    this->notify(INIT, state);
    stats_ = Stats{};
    info_  = AlgoInfo{};

    auto& energyF    = energyFunction_;
    auto gradientF   = derivative(energyF);
    auto hessianF    = derivative(gradientF);
    decltype(auto) e = energyF(x);
    decltype(auto) g = gradientF(x);
    decltype(auto) h = hessianF(x);

    eta_.resizeLike(g);
    Heta_.resizeLike(g);
    truncatedConjugateGradient_.analyzePattern(h);
    stats_.energy   = e;
    stats_.gradNorm = norm(g);
    truncatedConjugateGradient_.analyzePattern(h);

    innerInfo_.Delta = settings_.Delta0;
    spdlog::info(
        "        | iter | inner_i |   rho |   energy | energy_p | energy_inc |  norm(g) |    Delta | norm(corr) | "
        "InnerBreakReason");
    spdlog::info("{:-^143}", "-");
    while (not stoppingCriterion(e)) {
      this->notify(ITERATION_STARTED, state);
      if (settings_.useRand) {
        if (stats_.outerIter == 0) {
          eta_.setRandom();
          while (eta_.dot(h * eta_) > innerInfo_.Delta * innerInfo_.Delta)
            eta_ *= eps_; // eps is sqrt(sqrt(maschine-precision))
        } else
          eta_.setZero();
      } else
        eta_.setZero();

      solveInnerProblem(g, h);
      stats_.innerIterSum += innerInfo_.numInnerIter;

      info_.stopReasonString = tcg_stop_reason_[static_cast<int>(innerInfo_.stop_tCG)];
      Heta_                  = h * eta_;
      if (settings_.useRand and stats_.outerIter == 0) {
        info_.used_cauchy            = false;
        info_.randomPredictionString = " Used Random correction predictor";
        info_.cauchystr              = "                  ";
        double tauC;
        // Check the curvature
        const Eigen::VectorXd Hg = h * g;
        const auto g_Hg          = g.dot(Hg);
        if (g_Hg <= 0)
          tauC = 1;
        else
          tauC = std::min(Dune::power(stats_.gradNorm, 3) / (innerInfo_.Delta * g_Hg), 1.0);

        // generate the Cauchy point.
        const Eigen::VectorXd etaC  = -tauC * innerInfo_.Delta / stats_.gradNorm * g;
        const Eigen::VectorXd HetaC = -tauC * innerInfo_.Delta / stats_.gradNorm * Hg;

        const double mdle  = stats_.energy + g.dot(eta_) + .5 * Heta_.dot(eta_);
        const double mdlec = stats_.energy + g.dot(etaC) + .5 * HetaC.dot(etaC);
        if (mdlec < mdle && stats_.outerIter == 0) {
          eta_              = etaC;
          Heta_             = HetaC;
          info_.used_cauchy = true;

          info_.cauchystr = " USED CAUCHY POINT!";
        }
      } else
        info_.cauchystr = "                  ";

      stats_.etaNorm = eta_.norm();

      updateFunction_(x, eta_);

      // Calculate energy of our proposed update step
      e                     = energyF(x);
      stats_.energyProposal = e;

      // Will we accept the proposal or not?
      // Check the performance of the quadratic model against the actual energy.
      auto rhonum = stats_.energy - stats_.energyProposal;
      auto rhoden = -eta_.dot(g + 0.5 * Heta_);

      /*  Close to convergence the proposed energy and the real energy almost coincide.
       *  Therefore, the performance check of our model becomes ill-conditioned
       *  The regularisation fixes this */
      const auto rhoReg = std::max(1.0, abs(stats_.energy)) * eps_ * settings_.rho_reg;
      rhonum            = rhonum + rhoReg;
      rhoden            = rhoden + rhoReg;

      const bool modelDecreased = rhoden > 0.0;

      if (!modelDecreased)
        info_.stopReasonString.append(", model did not decrease");

      stats_.rho = rhonum / rhoden;
      stats_.rho = stats_.rho < 0.0 ? -1.0 : stats_.rho; // move rho to the domain [-1.0,inf]

      info_.trstr = "   ";

      // measure if energy decreased
      const bool energyDecreased = Dune::FloatCmp::ge(stats_.energy - stats_.energyProposal, -1e-12);

      // If the model behaves badly or if the energy increased we reduce the trust region size
      if (stats_.rho < 1e-4 || not modelDecreased || std::isnan(stats_.rho) || not energyDecreased) {
        info_.trstr = "TR-";
        innerInfo_.Delta /= 4.0;
        info_.consecutive_TRplus = 0;
        info_.consecutive_TRminus++;
        if (info_.consecutive_TRminus >= 5 && settings_.verbosity >= 1) {
          info_.consecutive_TRminus = -std::numeric_limits<int>::infinity();
          spdlog::info(" +++ Detected many consecutive TR- (radius decreases).");
          spdlog::info(" +++ Consider decreasing options.Delta_bar by an order of magnitude.");
        }

      } else if (stats_.rho > 0.99 && (innerInfo_.stop_tCG == Eigen::TCGStopReason::negativeCurvature ||
                                       innerInfo_.stop_tCG == Eigen::TCGStopReason::exceededTrustRegion)) {
        info_.trstr               = "TR+";
        innerInfo_.Delta          = std::min(3.5 * innerInfo_.Delta, settings_.Delta_bar);
        info_.consecutive_TRminus = 0;
        info_.consecutive_TRplus++;
        if (info_.consecutive_TRplus >= 5 && settings_.verbosity >= 1) {
          info_.consecutive_TRplus = -std::numeric_limits<int>::infinity();
          spdlog::info(" +++ Detected many consecutive TR+ (radius increases)");
          spdlog::info(" +++ Consider increasing options.Delta_bar by an order of magnitude");

        } else {
          info_.consecutive_TRplus  = 0;
          info_.consecutive_TRminus = 0;
        }
      }

      if (modelDecreased && stats_.rho > settings_.rho_prime && energyDecreased) {
        if (stats_.energyProposal > stats_.energy)
          spdlog::info(
              "Energy function increased by {} (step size: {}). Since this is small we accept the step and hope for "
              "convergence of the gradient norm.",
              stats_.energyProposal - stats_.energy, stats_.etaNorm);

        info_.acceptProposal       = true;
        info_.accstr               = "acc";
        info_.consecutive_Rejected = 0;
      } else {
        info_.acceptProposal = false;
        info_.accstr         = "REJ";

        if (info_.consecutive_Rejected >= 5)
          innerInfo_.Delta /= 2;
        else
          innerInfo_.Delta = std::min(innerInfo_.Delta, stats_.etaNorm / 2.0);
        ++info_.consecutive_Rejected;
      }

      stats_.outerIter++;

      if (settings_.verbosity == 1)
        logState();

      info_.randomPredictionString = "";

      solverInformation.correctionNorm = stats_.etaNorm;
      solverInformation.residualNorm   = stats_.gradNorm;
      this->notify(CORRECTION_UPDATED, state);

      if (info_.acceptProposal) {
        stats_.energy = stats_.energyProposal;
        this->notify(CORRECTIONNORM_UPDATED, state);
        this->notify(RESIDUALNORM_UPDATED, state);
        this->notify(SOLUTION_CHANGED, state);
      } else {
        updateFunction_(x, -eta_);
        eta_.setZero();
      }
      e               = energyF(x);
      g               = gradientF(x);
      h               = hessianF(x);
      stats_.gradNorm = g.norm();
      this->notify(NonLinearSolverMessages::ITERATION_ENDED, state);
    }
    spdlog::info("{}", info_.reasonString);
    spdlog::info("Total iterations: {} Total CG Iterations: {}", stats_.outerIter, stats_.innerIterSum);

    solverInformation.success =
        (info_.stop == StopReason::correctionNormTolReached) or (info_.stop == StopReason::gradientNormTolReached);

    solverInformation.iterations   = stats_.outerIter;
    solverInformation.residualNorm = stats_.gradNorm;
    if (solverInformation.success)
      this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, state);
    return solverInformation;
  }

  /**
   * \brief Access the energy function.
   * \return Reference to the energy function.
   */
  auto& energy() { return energyFunction_; }

  /**
   * \brief Access the energy function.
   * \return Reference to the energy function.
   */
  const auto& energy() const { return energyFunction_; }

  /**
   * \brief Access the residual.
   * \return The residual by value.
   */
  auto residual() const { return derivative(energyFunction_); }

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
  template <class T>
  constexpr auto make_optional_reference(T& value) {
    return std::make_optional<std::reference_wrapper<const T>>(std::cref(value));
  }

  template <class T>
  requires(not std::is_lvalue_reference_v<T>)
  constexpr T make_optional_reference(T&& value) {
    return value;
  }

  void logState() const {
    spdlog::info(
        "{:>3s} {:>3s} {:>6d} {:>9d}  {:>6.2f}  {:>9.2e}  {:>9.2e}  {:>11.2e}  {:>9.2e}  {:>9.2e}  {:>11.2e}   "
        "{:<73}",
        info_.accstr, info_.trstr, stats_.outerIter, innerInfo_.numInnerIter, stats_.rho, stats_.energy,
        stats_.energyProposal, stats_.energyProposal - stats_.energy, stats_.gradNorm, innerInfo_.Delta, stats_.etaNorm,
        info_.stopReasonString + info_.cauchystr + info_.randomPredictionString);
  }

  void logFinalState() {
    spdlog::info("{:>3s} {:>3s} {:>6d} {:>9d}  {: ^6}  {: ^9}  {: ^9}  {: ^11}  {:>9.2e}  {: ^9}  {: ^11}   {:<73}",
                 info_.accstr, info_.trstr, stats_.outerIter, innerInfo_.numInnerIter, " ", " ", " ", " ",
                 stats_.gradNorm, " ", " ", info_.stopReasonString + info_.cauchystr + info_.randomPredictionString);
  }

  bool stoppingCriterion(const auto& energy) {
    std::ostringstream stream;
    /** Gradient correction tolerance reached  */
    if (stats_.gradNorm < settings_.grad_tol && stats_.outerIter != 0) {
      logFinalState();
      spdlog::info("CONVERGENCE:  Energy: {:1.16e}    norm(gradient): {:1.16e}", energy, stats_.gradNorm);
      stream << "Gradient norm tolerance reached; options.tolerance = " << settings_.grad_tol;

      info_.reasonString = stream.str();

      info_.stop = StopReason::gradientNormTolReached;
      return true;
    } else if (stats_.etaNorm < settings_.corr_tol && stats_.outerIter != 0) {
      logFinalState();
      spdlog::info("CONVERGENCE:  Energy: {:1.16e}    norm(correction): {:1.16e}", energy, stats_.etaNorm);
      stream << "Displacement norm tolerance reached;  = " << settings_.corr_tol << "." << std::endl;

      info_.reasonString = stream.str();
      info_.stop         = StopReason::correctionNormTolReached;
      return true;
    }

    /** Maximum Time reached  */
    if (stats_.time >= settings_.maxtime) {
      logFinalState();
      stream << "Max time exceeded; options.maxtime = " << settings_.maxtime << ".";
      info_.reasonString = stream.str();
      info_.stop         = StopReason::maximumTimeReached;
      return true;
    }

    /** Maximum Iterations reached  */
    if (stats_.outerIter >= settings_.maxIter) {
      logFinalState();
      stream << "Max iteration count reached; options.maxiter = " << settings_.maxIter << ".";
      info_.reasonString = stream.str();
      info_.stop         = StopReason::maximumIterationsReached;
      return true;
    }
    return false;
  }

  void solveInnerProblem(const auto& g, const auto& h) {
    truncatedConjugateGradient_.setInfo(innerInfo_);
    int attempts = 0;
    truncatedConjugateGradient_.factorize(h);
    // If the preconditioner is IncompleteCholesky the factorization may fail if we have negative diagonal entries and
    // the initial shift is too small. Therefore, if the factorization fails we increase the initial shift by a factor
    // of 5.
    if constexpr (preConditioner == PreConditioner::IncompleteCholesky) {
      while (truncatedConjugateGradient_.info() != Eigen::Success) {
        choleskyInitialShift_ *= 5;
        truncatedConjugateGradient_.preconditioner().setInitialShift(choleskyInitialShift_);
        truncatedConjugateGradient_.factorize(h);
        if (attempts > 5)
          DUNE_THROW(Dune::MathError, "Factorization of preconditioner failed!");
        ++attempts;
      }
      if (truncatedConjugateGradient_.info() == Eigen::Success)
        choleskyInitialShift_ = 1e-3;
    }
    eta_       = truncatedConjugateGradient_.solveWithGuess(-g, eta_);
    innerInfo_ = truncatedConjugateGradient_.getInfo();
  }

  F energyFunction_;

  UpdateFunction updateFunction_;
  IDBCForceFunction idbcForceFunction_;
  std::remove_cvref_t<CorrectionType> eta_;
  std::remove_cvref_t<CorrectionType> Heta_;
  Settings settings_;
  AlgoInfo info_;
  double choleskyInitialShift_ = 1e-3;
  Eigen::TCGInfo<double> innerInfo_;
  Stats stats_;
  static constexpr double eps_ = 0.0001220703125; // 0.0001220703125 is sqrt(sqrt(maschine-precision))
  std::array<std::string, 6> tcg_stop_reason_{
      {"negative curvature", "exceeded trust region", "reached target residual-kappa (linear)",
       "reached target residual-theta (superlinear)", "maximum inner iterations", "model increased"}
  };

  using PreConditionerType =
      std::conditional_t<preConditioner == PreConditioner::IdentityPreconditioner, Eigen::IdentityPreconditioner,
                         std::conditional_t<preConditioner == PreConditioner::DiagonalPreconditioner,
                                            typename Eigen::DiagonalPreconditioner<std::decay_t<EnergyType>>,
                                            typename Eigen::IncompleteCholesky<std::decay_t<EnergyType>>>>;
  Eigen::TruncatedConjugateGradient<std::decay_t<HessianType>, Eigen::Lower | Eigen::Upper, PreConditionerType>
      truncatedConjugateGradient_;
};

/**
 * \brief Creates an instance of the TrustRegion solver.
 *
 * \tparam F Type of the function to solve.
 * \tparam preConditioner Type of the preconditioner used internally (default is IncompleteCholesky).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \param f The function to solve.
 * \param updateFunction Update function (default is UpdateDefault).
 * \return Shared pointer to the TrustRegion solver instance.
 */
template <typename F, PreConditioner preConditioner = PreConditioner::IncompleteCholesky,
          typename UF = utils::UpdateDefault, typename IDBCF = utils::IDBCForceDefault>
auto makeTrustRegion(const F& f, UF&& updateFunction = {}, IDBCF&& idbcForceFunction = {}) {
  return std::make_shared<TrustRegion<F, preConditioner, UF, IDBCF>>(f, updateFunction, idbcForceFunction);
}

template <typename F, PreConditioner preConditioner = PreConditioner::IncompleteCholesky,
          typename UF2 = utils::UpdateDefault, typename IDBCF2 = utils::IDBCForceDefault>
TrustRegion(const F& f, UF2&& updateFunction = {}, IDBCF2&& idbcForceFunction = {})
    -> TrustRegion<F, preConditioner, std::remove_cvref_t<UF2>, std::remove_cvref_t<IDBCF2>>;

} // namespace Ikarus
