// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file functionsanitychecks.hh
 * \brief Implementation of function sanity checks
 */

#pragma once

#include <type_traits>

#include <dune/common/float_cmp.hh>
#include <dune/functions/common/signature.hh>

#include <spdlog/spdlog.h>

#include <ikarus/utils/defaultfunctions.hh>

namespace Ikarus::utils {
namespace Impl {
  /**
   * \internal
   * \brief Draws the function and returns the slope for a given function.
   * \param functionName The name of the function.
   * \param ftfunc The function to be evaluated.
   * \param draw Flag indicating whether to draw the result.
   * \param slopeOfReference The reference slope for comparison.
   * \return The computed slope.
   */
  double drawResultAndReturnSlope(std::string&& functionName, const std::function<double(double)>& ftfunc, bool draw,
                                  int slopeOfReference);
} // namespace Impl
/**
 * \brief Struct to hold flags for function checks.
 */
struct CheckFlags
{
  bool draw                        = true;
  bool writeSlopeStatementIfFailed = true;
  double tolerance                 = 1e-2;
};

/**
 * \brief Checks the gradient of a differentiable Functions.
 * \details The checkgradient function is inspired by http://sma.epfl.ch/~nboumal/book/  Chapter 4.8 and
 * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/checkdiff.m
 * \ingroup utils
 * \tparam F Type of the differentiable Functions.
 * \tparam UF Type of the update function.
 * \param f The differentiable Functions.
 * \param x the evaluation point.
 * \param checkFlags Flags for the check.
 * \param p_updateFunction Update function.
 * \return True if the check passed, false otherwise.
 */
template <typename F, typename UF = UpdateDefault>
bool checkGradient(F& f, const typename F::Domain& p, CheckFlags checkFlags = CheckFlags(),
                   UF&& p_updateFunction = {}) {
  auto x           = p;
  auto gradF       = derivative(f);
  const auto e     = f(x);
  decltype(auto) g = gradF(x);
  using UpdateType = std::remove_cvref_t<decltype(g)>;

  UpdateType b;
  if constexpr (not std::is_floating_point_v<UpdateType>) {
    b.resizeLike(g);
    b.setRandom();
    b /= b.norm();
  } else
    b = 1;

  double gradfv;
  if constexpr (not std::is_floating_point_v<UpdateType>)
    gradfv = g.dot(b);
  else
    gradfv = g * b;

  auto ftfunc = [&](auto t) {
    p_updateFunction(x, t * b);
    const auto ept = f(x);
    auto value     = std::abs(ept - e - t * gradfv);
    p_updateFunction(x, -t * b);
    return value;
  };

  const double slope = Impl::drawResultAndReturnSlope("Gradient", ftfunc, checkFlags.draw, 2);

  const bool checkPassed = Dune::FloatCmp::le(2.0, slope, checkFlags.tolerance);

  if (checkFlags.writeSlopeStatementIfFailed and not checkPassed) {
    spdlog::info("Gradient check:");
    spdlog::info("The slope should be 2. It seems to be {}.", slope);
    if (checkPassed)
      spdlog::info("We consider this as sufficient.");
    else
      spdlog::info("The gradient seems wrong.");
  }

  return checkPassed;
}

/**
 * \brief Checks the Jacobian of a differentiable Functions.
 * \details The checkjacobian function is inspired by http://sma.epfl.ch/~nboumal/book/  Chapter 4.8 and
 * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/checkdiff.m
 * \ingroup utils
 * \tparam F Type of the differentiable Functions.
 * \tparam UF Type of the update function.
 * \param f The differentiable Functions.
 * \param x the evaluation point.
 * \param checkFlags Flags for the check.
 * \param p_updateFunction Update function.
 * \return True if the check passed, false otherwise.
 */
template <typename F, typename UF = UpdateDefault>
bool checkJacobian(F& f, const typename F::Domain& p, CheckFlags checkFlags = CheckFlags(),
                   UF&& p_updateFunction = {}) {
  auto x = p;

  auto gradF       = derivative(f);
  const auto e     = f(x);
  decltype(auto) g = gradF(x);
  using UpdateType =
      std::conditional_t<F::nDerivatives == 2, std::remove_cvref_t<decltype(g)>, std::remove_cvref_t<decltype(e)>>;

  UpdateType b;
  b.resizeLike(g.col(0));
  b.setRandom();
  b /= b.norm();

  const auto jacofv = (g * b).eval();

  auto ftfunc = [&](auto t) {
    p_updateFunction(x, t * b);
    const auto etb = f(x);
    auto value     = (etb - e - t * jacofv).norm();
    p_updateFunction(x, -t * b);
    return value;
  };

  const double slope = Impl::drawResultAndReturnSlope("Jacobian", ftfunc, checkFlags.draw, 2);

  const bool checkPassed = Dune::FloatCmp::le(2.0, slope, checkFlags.tolerance);

  if (checkFlags.writeSlopeStatementIfFailed and not checkPassed) {
    spdlog::info("Jacobian check:");
    spdlog::info("The slope should be 2. It seems to be {}.", slope);
    if (checkPassed)
      spdlog::info("We consider this as sufficient.");
    else
      spdlog::info("The Jacobian seems wrong.");
  }
  return checkPassed;
}

/**
 * \brief Checks the Hessian of a differentiable Functions.
 * \details  The checkHessian function is inspired by http://sma.epfl.ch/~nboumal/book/  Chapter 6.8 and
 * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/checkhessian.m
 * \ingroup utils
 * \tparam F Type of the differentiable Functions.
 * \tparam UF Type of the update function.
 * \param f The differentiable Functions.
 * \param x the evaluation point.
 * \param checkFlags Flags for the check.
 * \param p_updateFunction Update function.
 * \return True if the check passed, false otherwise.
 */
template <typename F, typename UF = UpdateDefault>
bool checkHessian(F& f, const typename F::Domain& p, CheckFlags checkFlags = CheckFlags(), UF&& p_updateFunction = {}) {
  auto x           = p;
  auto gradF       = derivative(f);
  auto hessF       = derivative(gradF);
  const auto e     = f(x);
  decltype(auto) g = gradF(x);
  decltype(auto) h = hessF(x);
  using UpdateType = std::remove_cvref_t<decltype(g)>;
  UpdateType b;

  if constexpr (not std::is_floating_point_v<UpdateType>) {
    b.resizeLike(g);
    b.setRandom();
    b /= b.norm();
  } else
    b = 1;

  double gradfv, vhessv;
  if constexpr (not std::is_floating_point_v<UpdateType>) {
    gradfv = g.dot(b);
    vhessv = (h * b).dot(b);
  } else {
    gradfv = g * b;
    vhessv = h * b * b;
  }

  auto ftfunc = [&](auto t) {
    p_updateFunction(x, t * b);
    const auto etb = f(x);
    auto value     = std::abs(etb - e - t * gradfv - 0.5 * t * t * vhessv);
    p_updateFunction(x, -t * b);
    return value;
  };

  const double slope = Impl::drawResultAndReturnSlope("Hessian", ftfunc, checkFlags.draw, 3);

  const bool checkPassed = Dune::FloatCmp::le(3.0, slope, checkFlags.tolerance);

  if (checkFlags.writeSlopeStatementIfFailed and not checkPassed) {
    spdlog::info("Hessian check:");
    spdlog::info("The slope should be 3. It seems to be {}.", slope);
    if (checkPassed)
      spdlog::info("We consider this as sufficient.");
    else
      spdlog::info("The Hessian seems wrong.");
  }
  return checkPassed;
}
} // namespace Ikarus::utils
