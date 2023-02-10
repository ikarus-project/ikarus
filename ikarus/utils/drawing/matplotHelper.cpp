// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

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
