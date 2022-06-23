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
#include <ikarus/localFunctions/expressions/binaryExpr.hh>
#include <ikarus/localFunctions/expressions/rebind.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename E1, typename E2>
  class LocalFunctionSum : public BinaryLocalFunctionExpression<LocalFunctionSum, E1, E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionSum, E1, E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionSum>;

    template <size_t ID_ = 0>
    static constexpr int orderID = std::max(Base::E1Raw::template order<ID_>(), Base::E1Raw::template order<ID_>());

    using ctype                    = typename Traits::ctype;
    static constexpr int valueSize = Traits::valueSize;
    static constexpr int gridDim   = Traits::gridDim;

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return Ikarus::eval(evaluateFunctionImpl(this->l(), localFunctionArgs)
                          + evaluateFunctionImpl(this->r(), localFunctionArgs));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return Ikarus::eval(evaluateDerivativeImpl(this->l(), localFunctionArgs)
                          + evaluateDerivativeImpl(this->r(), localFunctionArgs));
    }
  };

  template <typename E1, typename E2>
  struct LocalFunctionTraits<LocalFunctionSum<E1, E2>> : public LocalFunctionTraits<std::remove_cvref_t<E1>> {
    using Base = LocalFunctionTraits<std::remove_cvref_t<E1>>;
  };

  template <typename E1, typename E2>
  requires IsLocalFunction<E1, E2>
  constexpr auto operator+(E1&& u, E2&& v) {
    static_assert(Concepts::AddAble<typename std::remove_cvref_t<E1>::FunctionReturnType,
                                    typename std::remove_cvref_t<E2>::FunctionReturnType>,
                  "The function values of your local functions are not addable!");
    return LocalFunctionSum<E1, E2>(std::forward<E1>(u), std::forward<E2>(v));
  }
}  // namespace Ikarus