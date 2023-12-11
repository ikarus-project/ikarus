// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>

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
