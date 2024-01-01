// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>
#include <iomanip>
#include <iostream>
#include <vector>

#include <dune/common/float_cmp.hh>

namespace Eigen {
  template <typename Derived>
  struct EigenBase;
}

template <typename Derived, typename OtherDerived>
requires(std::convertible_to<Derived, Eigen::EigenBase<Derived> const&>and std::convertible_to<
         OtherDerived, Eigen::EigenBase<OtherDerived> const&>) bool isApproxSame(Derived const& val,
                                                                                 OtherDerived const& other,
                                                                                 double prec) {
  if constexpr (requires {
                  val.isApprox(other, prec);
                  (val - other).isMuchSmallerThan(1, prec);
                })
    return val.isApprox(other, prec) or (val - other).isZero(prec);
  else if constexpr (requires { val.isApprox(other, prec); })
    return val.isApprox(other, prec);
  else  // Eigen::DiagonalMatrix branch
    return val.diagonal().isApprox(other.diagonal(), prec) or (val.diagonal() - other.diagonal()).isZero(prec);
}

template <typename TestSuiteType, typename ScalarType>
void checkScalars(TestSuiteType& t, const ScalarType val, const ScalarType expectedVal,
                  const std::string& messageIfFailed = "") {
  if constexpr (std::is_integral_v<ScalarType>)
    t.check(val == expectedVal) << std::setprecision(16) << "Incorrect Scalars:\t" << expectedVal << " Actual:\t" << val
                                << messageIfFailed;
  else
    t.check(Dune::FloatCmp::eq(val, expectedVal))
        << std::setprecision(16) << "Incorrect Scalars:\t" << expectedVal << " Actual:\t" << val << messageIfFailed;
}

template <typename TestSuiteType, typename ControlInformation>
void checkSolverInfos(TestSuiteType& t, const std::vector<int>& expectedIterations,
                      const ControlInformation& controlInfo, const int loadSteps,
                      const std::string& messageIfFailed = "") {
  for (size_t i = 0U; i < loadSteps; ++i) {
    t.check(expectedIterations[i] == controlInfo.solverInfos[i].iterations)
        << "Incorrect number of iterations at step " << i << " with expected:\t" << expectedIterations[i]
        << " and actual:\t" << controlInfo.solverInfos[i].iterations << messageIfFailed;
    t.check(controlInfo.solverInfos[i].success) << "Failed to converge at step " << i;
  }
}
