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

#include <ikarus/localFunctions/expressions/scalarUnaryExpressions/scalarUnaryExpression.hh>

namespace Ikarus {

  struct SqrtFunc {
    template <typename ScalarType>
    static auto value(const ScalarType& v) {
      return sqrt(v);
    }

    template <typename ScalarType>
    static auto derivative(const ScalarType& v) {
      return ScalarType(1.0) / (ScalarType(2.0) * sqrt(v));
    }

    template <typename ScalarType>
    static auto secondDerivative(const ScalarType& v) {
      using std::pow;
      return ScalarType(-1.0) / (ScalarType(4.0) * pow(v, 3.0 / 2.0));
    }

    template <typename ScalarType>
    static auto thirdDerivative(const ScalarType& v) {
      using std::pow;
      return ScalarType(3.0) / (ScalarType(8.0) * pow(v, 5.0 / 2.0));
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
