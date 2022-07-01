/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#pragma once
#include <matplot/matplot.h>

#include <Eigen/Core>
namespace Ikarus::plot {
  void draw_xy(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

  template <typename FunctionType>
  void drawFunction(FunctionType&& f, std::pair<double, double>&& xRange, int eValuationPoints = 100) {
    std::vector<double> x = matplot::linspace(xRange.first, xRange.second, eValuationPoints);
    std::vector<double> y = matplot::transform(x, [&f](auto x) { return f(x); });
    matplot::plot(x, y, "-o");
    matplot::hold(matplot::on);
  }

}  // namespace Ikarus::plot