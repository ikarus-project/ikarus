//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
namespace Ikarus {

  template <typename E1>
  class LocalFunctionNegate: public UnaryLocalFunctionExpression<LocalFunctionNegate<E1>, E1> {
  public:
    using Base = UnaryLocalFunctionExpression<LocalFunctionNegate<E1>, E1>;
    using Base::UnaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionNegate>;

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return eval(-evaluateFunctionImpl(this->_u, localFunctionArgs));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return eval(-evaluateDerivativeImpl(this->_u, localFunctionArgs));
    }
  };

  template <typename E1>
  struct LocalFunctionTraits<LocalFunctionNegate<E1>> {
    static constexpr int valueSize = E1::valueSize;
    using DomainType = typename E1::DomainType;
  };

  template <typename E1>
  LocalFunctionNegate<E1> operator-(LocalFunctionInterface<E1> const& u) {
    return LocalFunctionNegate<E1>(u);
  }

}  // namespace Ikarus