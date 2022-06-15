//
// Created by lex on 6/14/22.
//

#pragma once

#include <ikarus/localFunctions/expressions/scalarUnaryExpressions/scalarUnaryExpression.hh>

namespace Ikarus {

  struct LogFunc {
    template <typename ScalarType>
    static auto value(const ScalarType& v) {
      return log(v);
    }

    template <typename ScalarType>
    static auto derivative(const ScalarType& v) {
      return ScalarType(1) / v;
    }

    template <typename ScalarType>
    static auto secondDerivative(const ScalarType& v) {
      return ScalarType(-1) / (v * v);
    }

    template <typename ScalarType>
    static auto thirdDerivative(const ScalarType& v) {
      return ScalarType(2) / (v * v * v);
    }
  };
  template <typename E1>
  requires IsLocalFunction<E1>
  constexpr auto log(E1&& u) {
    static_assert(std::remove_cvref_t<E1>::valueSize == 1,
                  "Log expression only defined for scalar valued local functions.");
    return ScalarUnaryExpression<E1, LogFunc>(std::forward<E1>(u));
  }
}  // namespace Ikarus
