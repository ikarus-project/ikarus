//
// Created by alex on 3/8/22.
//


#pragma once
#include <functional>
#include <matplot/matplot.h>
#include <matplot/util/colors.h>

#include <dune/common/float_cmp.hh>

#include <ikarus/utils/utils/findLineSegment.h>

/*
 * The checkgradient function is inspired by http://sma.epfl.ch/~nboumal/book/  Chapter 4.8 and
 * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/checkdiff.m
 */
template <typename NonlinearOperator, typename UpdateType = typename NonlinearOperator::template Parameter<0>>
bool checkGradient(
    NonlinearOperator& nonLinOp, bool draw = true,
    std::function<void(typename NonlinearOperator::template Parameter<0>&, const UpdateType&)> p_updateFunction
    = [](typename NonlinearOperator::template Parameter<0>& a, const UpdateType& b) { a += b; }) {
  auto& x         = nonLinOp.firstParameter();
  const auto xOld = x;
  UpdateType b;
  if constexpr (not std::is_floating_point_v<UpdateType>) {
    b.resizeLike(nonLinOp.derivative());
    b.setRandom();
    b /= b.norm();
  } else
    b = 1;

  nonLinOp.template updateAll();
  const auto e = nonLinOp.value();

  double gradfv;
  if constexpr (not std::is_floating_point_v<UpdateType>)
    gradfv = nonLinOp.derivative().dot(b);
  else
    gradfv = nonLinOp.derivative() * b;

  auto ftfunc = [&](auto t) {
    p_updateFunction(x, t * b);
    nonLinOp.template update<0>();
    auto value = std::abs(nonLinOp.value() - e - t * gradfv);
    x          = xOld;
    return value;
  };

  using namespace matplot;
  auto f   = figure(true);
  auto ax1 = gca();
  hold(ax1, true);
  std::vector<double> t = logspace(-8, 0, 100);
  Eigen::Map<Eigen::VectorXd> tE(t.data(), t.size());
  std::vector<double> ftevaluated = transform(t, ftfunc);
  Eigen::Map<Eigen::VectorXd> yE(ftevaluated.data(), ftevaluated.size());

  std::vector<double> fexpectedSlope = transform(t, [](auto t) { return t * t; });
  const int rangeSize                = 10;
  const auto [poly, range]           = Ikarus::findLineSegment(tE.array().log10(), yE.array().log10(), rangeSize);
  spdlog::info("Gradient check:");
  spdlog::info("The slope should be 2. It seems to be {}.", poly.coefficients()[1]);
  bool checkPassed = Dune::FloatCmp::eq(2.0, poly.coefficients()[1], 1e-4);
  if (checkPassed)
    spdlog::info("We consider this as sufficient.");
  else
    spdlog::info("The gradient seems wrong.");

  if (draw) {
    std::vector<double> tOfRange(rangeSize);
    std::vector<double> fInRange(rangeSize);
    auto tET = tE(range);
    auto yET = yE(range);

    for (int i = 0; auto r : tET) {
      tOfRange[i] = r;
      fInRange[i] = yET[i];
      ++i;
    }

    auto l0          = ax1->loglog(t, ftevaluated);
    auto lexpected   = ax1->loglog(t, fexpectedSlope, "--");
    auto lFoundRange = ax1->loglog(tOfRange, fInRange);
    l0->line_width(2);
    lexpected->line_width(2);
    lFoundRange->line_width(4);
    lFoundRange->color("magenta");
    l0->color("blue");
    lexpected->color("red");
    xlabel("h");
    ylabel("Approximation error ");
    title("Gradient check");
    f->show();
  }
  return checkPassed;
}

/*
 * The checkgradient function is inspired by http://sma.epfl.ch/~nboumal/book/  Chapter 6.8 and
 * https://github.com/NicolasBoumal/manopt/blob/master/manopt/tools/checkhessian.m
 */
template <typename NonlinearOperator, typename UpdateType = typename NonlinearOperator::template Parameter<0>>
bool checkHessian(
    NonlinearOperator& nonLinOp, bool draw = true,
    std::function<void(typename NonlinearOperator::template Parameter<0>&, const UpdateType&)> p_updateFunction
    = [](typename NonlinearOperator::template Parameter<0>& a, const UpdateType& b) { a += b; }) {
  auto& x         = nonLinOp.firstParameter();
  const auto xOld = x;
  UpdateType b;
  if constexpr (not std::is_floating_point_v<UpdateType>) {
    b.resizeLike(nonLinOp.derivative());
    b.setRandom();
    b /= b.norm();
  } else
    b = 1;

  nonLinOp.template updateAll();
  const auto e = nonLinOp.value();

  double gradfv, vhessv;
  if constexpr (not std::is_floating_point_v<UpdateType>) {
    gradfv = nonLinOp.derivative().dot(b);
    vhessv = (nonLinOp.secondDerivative() * b).dot(b);
  } else {
    gradfv = nonLinOp.derivative() * b;
    vhessv = nonLinOp.secondDerivative() * b * b;
  }

  auto ftfunc = [&](auto t) {
    p_updateFunction(x, t * b);
    nonLinOp.template update<0>();
    auto value = std::abs(nonLinOp.value() - e - t * gradfv - 0.5 * t * t * vhessv);
    x          = xOld;
    return value;
  };
  using namespace matplot;
  auto f   = figure(true);
  auto ax1 = gca();
  hold(ax1, true);
  std::vector<double> t = logspace(-8, 0, 100);
  Eigen::Map<Eigen::VectorXd> tE(t.data(), t.size());
  std::vector<double> ftevaluated = transform(t, ftfunc);
  Eigen::Map<Eigen::VectorXd> yE(ftevaluated.data(), ftevaluated.size());

  std::vector<double> fexpectedSlope = transform(t, [](auto t) { return t * t * t; });
  const int rangeSize                = 10;
  const auto [poly, range]           = Ikarus::findLineSegment(tE.array().log10(), yE.array().log10(), rangeSize);
  spdlog::info("Hessian check:");
  spdlog::info("The slope should be 3. It seems to be {}.", poly.coefficients()[1]);
  bool checkPassed = Dune::FloatCmp::eq(3.0, poly.coefficients()[1], 1e-3);

  if (checkPassed)
    spdlog::info("We consider this as sufficient.");
  else
    spdlog::info("The Hessian seems wrong.");

  if (draw) {
    std::vector<double> tOfRange(rangeSize);
    std::vector<double> fInRange(rangeSize);
    auto tET = tE(range);
    auto yET = yE(range);

    for (int i = 0; auto r : tET) {
      tOfRange[i] = r;
      fInRange[i] = yET[i];
      ++i;
    }

    auto l0          = ax1->loglog(t, ftevaluated);
    auto lexpected   = ax1->loglog(t, fexpectedSlope, "--");
    auto lFoundRange = ax1->loglog(tOfRange, fInRange);
    l0->line_width(2);
    lexpected->line_width(2);
    lFoundRange->line_width(4);
    lFoundRange->color("magenta");
    l0->color("blue");
    lexpected->color("red");
    xlabel("h");
    ylabel("Approximation error ");
    title("Hessian check");
    f->show();
  }
  return checkPassed;
}