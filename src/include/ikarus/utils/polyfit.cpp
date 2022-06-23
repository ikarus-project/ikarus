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


#include "polyfit.hh"
namespace Ikarus {

  /*
   * This function returns the polynom fitted onto the data pased in in the least square sense.
   * It also return the least square error.
   */
  std::tuple<Dune::Functions::Polynomial<double>, double> polyfit(const Eigen::Ref<const Eigen::VectorXd>& x,
                                                                  const Eigen::Ref<const Eigen::VectorXd>& y,
                                                                  const int deg) {
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(x.size(), deg + 1);
    for (int j = 1; j < deg + 1; ++j)
      A.col(j) = A.col(j - 1).cwiseProduct(x);

    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(A.rows(), A.cols());
    qr.compute(A);
    Eigen::VectorXd coeffs = qr.solve(y);

    std::vector<double> coeffsSTD(coeffs.begin(), coeffs.end());
    Dune::Functions::Polynomial<double> poly(std::move(coeffsSTD));
    return std::make_tuple(poly, (A * coeffs - y).norm());
  }
}  // namespace Ikarus