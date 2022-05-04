//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
namespace Ikarus {

  template <typename E1>
  class LocalFunctionSqrt: public UnaryLocalFunctionExpression<LocalFunctionSqrt, E1> {
  public:
    using Base = UnaryLocalFunctionExpression<LocalFunctionSqrt, E1>;
    using Base::UnaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionSqrt>;
    static constexpr int valueSize =  1;

    template<size_t ID_=0>
    static constexpr int order = nonLinear;

    template <typename LFArgs>
    auto evaluateValueOfExpression(const LFArgs& lfArgs) const {
      return Ikarus::eval(sqrt(evaluateFunctionImpl(this->m(), lfArgs).getValue()));
    }

    template <int DerivativeOrder, typename LFArgs>
    auto evaluateDerivativeOfExpression(const LFArgs& lfArgs) const {
      const auto u = evaluateFunctionImpl(this->m(), lfArgs).getValue();
      if constexpr (DerivativeOrder == 1)  // d(sqrt(u(x)))/(dxdy) =  u_x /(2*sqrt(u(x))
      {
        const auto u_x = evaluateDerivativeImpl(this->l(), lfArgs);
        return Ikarus::eval(u_x/(2*sqrt(u)));
      }else if constexpr (DerivativeOrder == 2)
      {
        const auto& [u_x,u_y] = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs);
      }

    }
  };

  template <typename E1>
  struct LocalFunctionTraits<LocalFunctionSqrt<E1>> : public LocalFunctionTraits<std::remove_cvref_t<E1>>  {
  };

  template <typename E1> requires IsLocalFunction<E1>
  constexpr auto sqrt(E1 && u) {
    static_assert(E1::valueSize==1,"Sqrt expression only defined for scalar valued local functions.")
    return LocalFunctionSqrt<E1>(std::forward<E1>(u));
  }

}  // namespace Ikarus