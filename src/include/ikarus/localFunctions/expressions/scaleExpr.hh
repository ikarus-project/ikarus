//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/binaryExpr.hh>
namespace Ikarus {

  template <typename E1,typename E2>
  class LocalFunctionScale: public BinaryLocalFunctionExpression<LocalFunctionScale<E1,E2>, E1,E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionScale<E1,E2>, E1,E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionScale>;

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return eval(this->l*evaluateFunctionImpl(this->r, localFunctionArgs));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return eval(this->l*evaluateDerivativeImpl(this->r, localFunctionArgs));
    }
  };

  template <typename E1,typename E2>
  struct LocalFunctionTraits<LocalFunctionScale<E1,E2>> {
    static constexpr int valueSize = E2::valueSize;
    using DomainType = typename E2::DomainType;
  };

  template <typename E1,typename E2> requires std::is_arithmetic_v<E1>
  constexpr LocalFunctionScale<E1,E2> operator*(const E1& factor, LocalFunctionInterface<E2> const& u) {
    return LocalFunctionScale<E1,E2>(factor,u);
  }

  template <typename E1,typename E2> requires std::is_arithmetic_v<E1>
  constexpr LocalFunctionScale<E1,E2> operator*(LocalFunctionInterface<E2> const& u,const E1& factor) {
    return factor*u;
  }

}  // namespace Ikarus