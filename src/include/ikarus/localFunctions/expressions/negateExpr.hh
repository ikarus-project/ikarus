//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
namespace Ikarus {

  template <typename E1>
  class LocalFunctionNegate: public UnaryLocalFunctionExpression<LocalFunctionNegate, E1> {
  public:
    using Base = UnaryLocalFunctionExpression<LocalFunctionNegate, E1>;
    using Base::UnaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionNegate>;
    static constexpr int valueSize =  Traits::valueSize;

    template<size_t ID_=0>
    static constexpr int order = std::remove_cvref_t<E1>::template order<ID_>;

    template <typename LFArgs>
    auto evaluateValueOfExpression(const LFArgs& lfArgs) const {
      return Ikarus::eval(-evaluateFunctionImpl(this->m(), lfArgs));
    }

    template <int DerivativeOrder, typename LFArgs>
    auto evaluateDerivativeOfExpression(const LFArgs& lfArgs) const {
      return Ikarus::eval(- evaluateDerivativeImpl(this->m(), lfArgs));
    }
  };

  template <typename E1>
  struct LocalFunctionTraits<LocalFunctionNegate<E1>> : public LocalFunctionTraits<std::remove_cvref_t<E1>>  {
  };

  template <typename E1> requires IsLocalFunction<E1>
  constexpr auto operator-(E1 && u) {
    return LocalFunctionNegate<E1>(std::forward<E1>(u));
  }

}  // namespace Ikarus