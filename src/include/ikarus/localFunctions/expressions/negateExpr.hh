// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
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
