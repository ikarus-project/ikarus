//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/binaryExpr.hh>
#include <ikarus/localFunctions/expressions/arithmeticExpr.hh>
namespace Ikarus {

  template <typename E1,typename E2>
  class LocalFunctionScale: public BinaryLocalFunctionExpression<LocalFunctionScale<E1,E2>, E1,E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionScale<E1,E2>, E1,E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionScale>;

    template<typename OtherType,typename IDTuple>
    using Rebind = Rebind<LocalFunctionScale,E1,E2,OtherType,IDTuple>;

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return eval(this->l().value()*evaluateFunctionImpl(this->r(), localFunctionArgs));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return eval(this->l().value()*evaluateDerivativeImpl(this->r(), localFunctionArgs));
    }
  };

  template <typename E1,typename E2>
  struct LocalFunctionTraits<LocalFunctionScale<E1,E2>> {
    static constexpr int valueSize = E2::valueSize;
    using DomainType = typename E2::DomainType;
  };

  template <typename E1,typename E2> requires (std::is_arithmetic_v<E1> and IsLocalFunction<E2>and ! IsScaleExpr<E2>)
  constexpr LocalFunctionScale<ArithmeticExpr<E1>,E2> operator*(const E1& factor, E2 const& u) {
    return LocalFunctionScale<ArithmeticExpr<E1>,E2>(ArithmeticExpr(factor),u);
  }

  template <typename E1,typename E2> requires (std::is_arithmetic_v<E1> and IsLocalFunction<E2> and ! IsScaleExpr<E2>)
  constexpr LocalFunctionScale<ArithmeticExpr<E1>,E2> operator*(E2 const& u,const E1& factor) {
    return factor*u;
  }

// Simplification if nested scale expression occur
template <typename E1,typename E2> requires (std::is_arithmetic_v<E1> and IsScaleExpr<E2>)
constexpr auto operator*(const E1& factor, const E2& u) {
  return LocalFunctionScale<ArithmeticExpr<E1>,std::remove_cvref_t<decltype(u.r())>>(ArithmeticExpr(factor*u.l().value()),u.r());
}

template <typename E1,typename E2> requires (std::is_arithmetic_v<E1> and IsScaleExpr<E2> and IsLocalFunction<E2>)
constexpr auto operator*(E2 const& u,const E1& factor) {
  return factor*u;
}

template <typename E1,typename E2> requires (std::is_arithmetic_v<E1> and IsLocalFunction<E2> and ! IsScaleExpr<E2>)
constexpr LocalFunctionScale<ArithmeticExpr<E1>,E2> operator/(E2 const& u,const E1& factor) {
  static_assert(std::floating_point<E1>,"Operator/ should only called with floating point types.");

  return (1/factor)*u;
}

}  // namespace Ikarus