//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/binaryExpr.hh>
#include <ikarus/localFunctions/expressions/constant.hh>
namespace Ikarus {

  template <typename E1, typename E2>
  class LocalFunctionScale : public BinaryLocalFunctionExpression<LocalFunctionScale, E1, E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionScale, E1, E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits                   = LocalFunctionTraits<LocalFunctionScale>;
    static constexpr int valueSize = Traits::valueSize;

    template <size_t ID_ = 0>
    static constexpr int orderID
        = std::min(Base::E1Raw::template orderID<ID_> + Base::E2Raw::template orderID<ID_>, nonLinear);

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return Ikarus::eval(this->l().value() * evaluateFunctionImpl(this->r(), localFunctionArgs));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return Ikarus::eval(this->l().value() * evaluateDerivativeImpl(this->r(), localFunctionArgs));
    }
  };

  template <typename E1, typename E2>
  struct LocalFunctionTraits<LocalFunctionScale<E1, E2>> : public LocalFunctionTraits<std::remove_cvref_t<E2>> {};

  template <typename E1, typename E2>
  requires(std::is_arithmetic_v<std::remove_cvref_t<E1>>and
               IsLocalFunction<E2> and !IsScaleExpr<E2>) constexpr LocalFunctionScale<ConstantExpr<E1>, E2>
  operator*(E1&& factor, E2&& u) {
    return LocalFunctionScale<ConstantExpr<E1>, E2>(ConstantExpr(factor), std::forward<E2>(u));
  }

  template <typename E1, typename E2>
  requires(std::is_arithmetic_v<std::remove_cvref_t<E1>>and
               IsLocalFunction<E2> and !IsScaleExpr<E2>) constexpr LocalFunctionScale<ConstantExpr<E1>, E2>
  operator*(E2&& u, E1&& factor) {
    return factor * u;
  }

  // Simplification if nested scale expression occur
  template <typename E1, typename E2>
  requires(std::is_arithmetic_v<std::remove_cvref_t<E1>>and IsScaleExpr<E2>) constexpr auto operator*(E1&& factor,
                                                                                                      E2&& u) {
    u.l().value() *= factor;
    return u;
  }

  template <typename E1, typename E2>
  requires(std::is_arithmetic_v<std::remove_cvref_t<E1>>and IsScaleExpr<E2>and IsLocalFunction<E2>) constexpr auto
  operator*(E2&& u, E1&& factor) {
    return operator*(std::forward<E1>(factor), std::forward<E2>(u));
  }

  // Division operator
  template <typename E1, typename E2>
  requires(std::is_arithmetic_v<std::remove_cvref_t<E1>>and IsLocalFunction<E2> and !IsScaleExpr<E2>) constexpr auto
  operator/(E2&& u, E1&& factor) {
    static_assert(std::floating_point<std::remove_cvref_t<E1>>,
                  "Operator/ should only called with floating point types.");

    return operator*(std::forward<E1>(1 / factor), std::forward<E2>(u));
  }

}  // namespace Ikarus