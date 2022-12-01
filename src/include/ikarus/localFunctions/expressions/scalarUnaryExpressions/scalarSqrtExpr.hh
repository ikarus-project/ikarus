// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <ikarus/localFunctions/expressions/scalarUnaryExpressions/scalarUnaryExpression.hh>

namespace Ikarus {

  struct SqrtFunc {
    template <typename ScalarType>
    static auto value(const ScalarType& v) {
      return sqrt(v);
    }

    template <typename ScalarType>
    static auto derivative(const ScalarType& v) {
      return ScalarType(0.5) / sqrt(v);
    }

    template <typename ScalarType>
    static auto secondDerivative(const ScalarType& v) {
      using std::pow;
      return ScalarType(-0.25) / pow(v, 1.5);
    }

    template <typename ScalarType>
    static auto thirdDerivative(const ScalarType& v) {
      using std::pow;
      return ScalarType(0.375) / pow(v, 2.5);
    }
  };
  template <typename E1>
  requires IsLocalFunction<E1>
  constexpr auto sqrt(E1&& u) {
    static_assert(std::remove_cvref_t<E1>::valueSize == 1,
                  "Sqrt expression only defined for scalar valued local functions.");
    return ScalarUnaryExpression<E1, SqrtFunc>(std::forward<E1>(u));
  }
}  // namespace Ikarus
