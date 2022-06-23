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



#include "matplotHelper.hh"

namespace Ikarus::plot {
  void draw_xy(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    assert(x.size() == y.size() && "The passed x and y vectors have to have the same size!");
    std::vector<double> xstd(x.begin(), x.end());
    std::vector<double> ystd(y.begin(), y.end());
    matplot::plot(x, y, "-o");
    matplot::hold(matplot::on);
  }

}  // namespace Ikarus::plot