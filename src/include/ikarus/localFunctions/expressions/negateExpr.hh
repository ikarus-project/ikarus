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
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
namespace Ikarus {

  template <typename E1>
  class LocalFunctionNegate : public UnaryLocalFunctionExpression<LocalFunctionNegate, E1> {
  public:
    using Base = UnaryLocalFunctionExpression<LocalFunctionNegate, E1>;
    using Base::UnaryLocalFunctionExpression;
    using Traits                   = LocalFunctionTraits<LocalFunctionNegate>;
    static constexpr int valueSize = Traits::valueSize;

    template <size_t ID_ = 0>
    static constexpr int orderID = Base::E1Raw::template orderID<ID_>;

    template <typename LFArgs>
    auto evaluateValueOfExpression(const LFArgs& lfArgs) const {
      return Ikarus::eval(-evaluateFunctionImpl(this->m(), lfArgs));
    }

    template <int DerivativeOrder, typename LFArgs>
    auto evaluateDerivativeOfExpression(const LFArgs& lfArgs) const {
      return Ikarus::eval(-evaluateDerivativeImpl(this->m(), lfArgs));
    }
  };

  template <typename E1>
  struct LocalFunctionTraits<LocalFunctionNegate<E1>> : public LocalFunctionTraits<std::remove_cvref_t<E1>> {};

  template <typename E1>
  requires IsLocalFunction<E1>
  constexpr auto operator-(E1&& u) { return LocalFunctionNegate<E1>(std::forward<E1>(u)); }

}  // namespace Ikarus