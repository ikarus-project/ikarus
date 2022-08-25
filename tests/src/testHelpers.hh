
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

#include <Eigen/Core>
#include <catch2/matchers/catch_matchers_templated.hpp>

template<typename Derived>
struct EigenApproxEqual : Catch::Matchers::MatcherGenericBase {

  EigenApproxEqual(Eigen::EigenBase<Derived> const& vMB, double prec):
                                           v{ vMB.derived() }, prec{prec}
  {}

  template<typename OtherDerived>
  bool match(Eigen::EigenBase<OtherDerived> const& otherMB) const {
    const OtherDerived& other = otherMB.derived();
    if constexpr (requires {
                    v.isApprox(other, prec);
                    (v - other).isMuchSmallerThan(1, prec);
                  })
      return v.isApprox(other, prec) or (v - other).isZero(prec);
    else if constexpr (requires { v.isApprox(other, prec); })
      return v.isApprox(other, prec);
    else  // Eigen::DiagonalMatrix branch
      return v.diagonal().isApprox(other.diagonal(), prec) or (v.diagonal() - other.diagonal()).isZero(prec);
  }

  std::string describe() const override {
    std::stringstream stringstream;

    if constexpr (std::convertible_to<Derived, const Eigen::MatrixBase<Derived>&>) {
      stringstream << v;
    } else {  // branch for Eigen::DiagonalMatrix
      using Scalar = typename Derived::Scalar;
      using namespace Eigen;
      constexpr int diag_size = EIGEN_SIZE_MIN_PREFER_DYNAMIC(Derived::RowsAtCompileTime, Derived::ColsAtCompileTime);
      stringstream << Eigen::Matrix<Scalar, diag_size, diag_size>(v.derived().diagonal().asDiagonal());
    }
    return "Equals: " + stringstream.str();
  }
private:
  Derived const& v;
  double prec;
};



//MATCHER_P2(EigenApproxEqual, expect, prec,
//          std::string(negation ? "isn't" : "is") + " approx equal to\n" + ::testing::PrintToString(expect)
//              + "\nwith precision " + ::testing::PrintToString(prec)) {
//
//}
//
//MATCHER_P(EigenExactEqual, expect,
//         std::string(negation ? "isn't" : "is") + " equal to" + ::testing::PrintToString(expect)) {
// return ((arg == expect) == true).all();
//}
//
