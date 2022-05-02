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

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return Ikarus::eval(-evaluateFunctionImpl(this->m(), localFunctionArgs));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return Ikarus::eval(- evaluateDerivativeImpl(this->m(), localFunctionArgs) );
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