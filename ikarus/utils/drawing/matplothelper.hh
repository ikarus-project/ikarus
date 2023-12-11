// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <matplot/matplot.h>

#include <Eigen/Core>
namespace Ikarus::plot {
  void draw_xy(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

  template <typename FunctionType>
  void drawFunction(FunctionType&& f, std::pair<double, double>&& xRange, int eValuationPoints = 100) {
    std::vector<double> x = matplot::linspace(xRange.first, xRange.second, eValuationPoints);
    std::vector<double> y = matplot::transform(x, [&f](auto x_) { return f(x_); });
    matplot::plot(x, y, "-o");
    matplot::hold(matplot::on);
  }

}  // namespace Ikarus::plot
