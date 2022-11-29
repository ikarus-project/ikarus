// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <ikarus/localFunctions/expressions/scalarUnaryExpressions/scalarUnaryExpression.hh>

namespace Ikarus {

  template <int exponent>
  struct PowFunc {
    template <typename ScalarType>
    static auto value(const ScalarType& v) {
      return Dune::power(v, exponent);
    }

    template <typename ScalarType>
    static auto derivative(const ScalarType& v) {
      return Dune::power(v, exponent - 1) * exponent;
    }

    template <typename ScalarType>
    static auto secondDerivative(const ScalarType& v) {
      return Dune::power(v, exponent - 2) * exponent * (exponent - 1);
    }

    template <typename ScalarType>
    static auto thirdDerivative(const ScalarType& v) {
      return Dune::power(v, exponent - 3) * exponent * (exponent - 1) * (exponent - 2);
    }
  };
  template <int exponent, typename E1>
  requires IsLocalFunction<E1>
  constexpr auto pow(E1&& u) {
    static_assert(std::remove_cvref_t<E1>::valueSize == 1,
                  "Pow expression only defined for scalar valued local functions.");
    return ScalarUnaryExpression<E1, PowFunc<exponent>>(std::forward<E1>(u));
  }
}  // namespace Ikarus
