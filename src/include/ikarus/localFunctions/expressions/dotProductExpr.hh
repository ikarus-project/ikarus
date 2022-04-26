//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/binaryExpr.hh>
namespace Ikarus {

  template <typename E1, typename E2>
  class LocalFunctionDot : public BinaryLocalFunctionExpression<LocalFunctionDot<E1, E2>, E1, E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionDot<E1, E2>, E1, E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionDot>;

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return evaluateFunctionImpl(this->_u, localFunctionArgs)
          .getValue()
          .dot(evaluateFunctionImpl(this->_v, localFunctionArgs).getValue());
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      static_assert(DerivativeOrder == 1 or DerivativeOrder == 2,
                    "For testing only first/second order derivatives are supported.");
      if constexpr (DerivativeOrder == 1)  // d(dot(u,v))/dx =  u_x * v+ u*v_x
        return eval((transpose(evaluateDerivativeImpl(this->_u, localFunctionArgs))
                     * evaluateFunctionImpl(this->_v, localFunctionArgs).getValue())
                        .transpose()
                    + evaluateFunctionImpl(this->_u, localFunctionArgs).getValue().transpose()
                          * evaluateDerivativeImpl(this->_v, localFunctionArgs));
      else if constexpr (DerivativeOrder
                         == 2) {  // dd(dot(u,v))/(dxdy) =  u_{x,y} * v + u_x*v_y + u_y* v_x + u * v_{x,y}
        const auto u = evaluateFunctionImpl(this->_u, localFunctionArgs).getValue();
        const auto v = evaluateFunctionImpl(this->_v, localFunctionArgs).getValue();

        const auto firstWRTDerivativeArgs  = localFunctionArgs.extractSpatialOrFirstWrtArg();
        const auto secondWRTDerivativeArgs = localFunctionArgs.extractSecondWrtArgOrFirstNonSpatial();

        if constexpr (LocalFunctionEvaluationArgs_::hasOneSpatialSingle
                      and LocalFunctionEvaluationArgs_::hasSingleCoeff) {
          const auto u_xy = evaluateDerivativeImpl(this->_u, localFunctionArgs);
          const auto v_xy = evaluateDerivativeImpl(this->_v, localFunctionArgs);
          const auto u_x  = evaluateDerivativeImpl(this->_u, firstWRTDerivativeArgs);
          const auto u_y  = evaluateDerivativeImpl(this->_u, secondWRTDerivativeArgs);
          const auto v_x  = evaluateDerivativeImpl(this->_v, firstWRTDerivativeArgs);
          const auto v_y  = evaluateDerivativeImpl(this->_v, secondWRTDerivativeArgs);

          return eval(transpose(v) * u_xy + transpose(u_x) * v_y + transpose(v_x) * u_y + transpose(u) * v_xy);
        } else if constexpr (LocalFunctionEvaluationArgs_::hasNoSpatial and LocalFunctionEvaluationArgs_::hasTwoCoeff) {
          const auto alonguArgs = LocalFunctionEvaluationArgs_::createWithAlong(localFunctionArgs, along(v));
          const auto alongvArgs = LocalFunctionEvaluationArgs_::createWithAlong(localFunctionArgs, along(u));

          const auto [coeff0Args, coeff1Args] = localFunctionArgs.extractWrtArgsTwoCoeffsToSingleCoeff();

          const auto u_xyAlongv = evaluateDerivativeImpl(this->_u, alongvArgs);
          const auto v_xyAlongu = evaluateDerivativeImpl(this->_v, alonguArgs);
          const auto u_x        = evaluateDerivativeImpl(this->_u, coeff0Args);
          const auto u_y        = evaluateDerivativeImpl(this->_u, coeff1Args);
          const auto v_x        = evaluateDerivativeImpl(this->_v, coeff0Args);
          const auto v_y        = evaluateDerivativeImpl(this->_v, coeff1Args);

          return eval(u_xyAlongv + transpose(u_x) * v_y + transpose(v_x) * u_y + v_xyAlongu);
        }
      }
    }
  };

  template <typename E1, typename E2>
  struct LocalFunctionTraits<LocalFunctionDot<E1, E2>> {
    /** \brief Size of the function value */
    static constexpr int valueSize = 1;
    /** \brief Type for the points for evaluation, usually the integration points */
    using DomainType = typename E1::DomainType;
  };

  template <typename E1, typename E2>
  LocalFunctionDot<E1, E2> dot(LocalFunctionInterface<E1> const& u, LocalFunctionInterface<E2> const& v) {
    return LocalFunctionDot<E1, E2>(u, v);
  }

}  // namespace Ikarus