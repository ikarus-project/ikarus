//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
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
        const auto u_x = evaluateDerivativeImpl(this->m(), lfArgs);
        return Ikarus::eval(u_x/(2*sqrt(u)));
      }else if constexpr (DerivativeOrder == 2)
      {
        const auto& [u_x,u_y] = evaluateFirstOrderDerivativesImpl(this->m(), lfArgs);
        const auto u_xy = evaluateDerivativeImpl(this->m(), lfArgs);

        return Ikarus::eval(-u_x*u_y/(4*std::pow(u,3/2)) + u_xy/(2*sqrt(u)));
      }else if constexpr (DerivativeOrder == 3)
      {
        const auto& [u_x,u_y,u_z] = evaluateFirstOrderDerivativesImpl(this->m(), lfArgs);
        const auto &[u_xy, u_xz] = evaluateSecondOrderDerivativesImpl(this->l(), lfArgs);

        const auto u_xyz = evaluateDerivativeImpl(this->m(), lfArgs);

        const auto uppow32 = std::pow(u,3.0/2.0);
        const auto argsForDyz = lfArgs.extractSecondWrtArgOrFirstNonSpatial();
        const auto u_yz = evaluateDerivativeImpl(this->m(), argsForDyz);

        return Ikarus::eval((3*u_x*u_y*u_z)/(8*std::pow(u,5.0/2.0)) - u_xz*u_y/(4*uppow32) - u_x*u_yz/(4*uppow32) - u_xy*u_z/(4*uppow32) + u_xyz/(2*sqrt(u)));
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