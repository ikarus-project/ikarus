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
      return evaluateFunctionImpl(this->l, localFunctionArgs)
          .getValue()
          .dot(evaluateFunctionImpl(this->r, localFunctionArgs).getValue());
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      static_assert(DerivativeOrder == 1 or DerivativeOrder == 2 or DerivativeOrder == 3,
                    "For testing only first/second order derivatives are supported.");
      if constexpr (DerivativeOrder == 1)  // d(dot(u,v))/dx =  u_x * v+ u*v_x
        return eval((transpose(evaluateDerivativeImpl(this->l, localFunctionArgs))
                     * evaluateFunctionImpl(this->r, localFunctionArgs).getValue())
                        .transpose()
                    + evaluateFunctionImpl(this->l, localFunctionArgs).getValue().transpose()
                          * evaluateDerivativeImpl(this->r, localFunctionArgs));
      else if constexpr (DerivativeOrder == 2 ) {  // dd(dot(u,v))/(dxdy) =  u_{x,y} * v + u_x*v_y + u_y* v_x + u * v_{x,y}
        const auto u = evaluateFunctionImpl(this->l, localFunctionArgs).getValue();
        const auto v = evaluateFunctionImpl(this->r, localFunctionArgs).getValue();

        const auto argsForDx  = localFunctionArgs.extractSpatialOrFirstWrtArg();
        const auto argsForDy = localFunctionArgs.extractSecondWrtArgOrFirstNonSpatial();

        if constexpr (LocalFunctionEvaluationArgs_::hasOneSpatialSingle
                      and LocalFunctionEvaluationArgs_::hasSingleCoeff) {
          const auto u_xy = evaluateDerivativeImpl(this->l, localFunctionArgs);
          const auto v_xy = evaluateDerivativeImpl(this->r, localFunctionArgs);
          const auto u_x  = evaluateDerivativeImpl(this->l, argsForDx);
          const auto u_y  = evaluateDerivativeImpl(this->l, argsForDy);
          const auto v_x  = evaluateDerivativeImpl(this->r, argsForDx);
          const auto v_y  = evaluateDerivativeImpl(this->r, argsForDy);

          return eval(transpose(v) * u_xy + transpose(u_x) * v_y + transpose(v_x) * u_y + transpose(u) * v_xy);
        } else if constexpr (LocalFunctionEvaluationArgs_::hasNoSpatial and LocalFunctionEvaluationArgs_::hasTwoCoeff) {
          const auto alonguArgs = createWithAlong(localFunctionArgs, along(v));
          const auto alongvArgs = createWithAlong(localFunctionArgs, along(u));

          const auto [coeff0Args, coeff1Args] = localFunctionArgs.extractWrtArgsTwoCoeffsToSingleCoeff();

          const auto u_xyAlongv = evaluateDerivativeImpl(this->l, alongvArgs);
          const auto v_xyAlongu = evaluateDerivativeImpl(this->r, alonguArgs);
          const auto u_x        = evaluateDerivativeImpl(this->l, coeff0Args);
          const auto u_y        = evaluateDerivativeImpl(this->l, coeff1Args);
          const auto v_x        = evaluateDerivativeImpl(this->r, coeff0Args);
          const auto v_y        = evaluateDerivativeImpl(this->r, coeff1Args);

          return eval(u_xyAlongv + transpose(u_x) * v_y + transpose(v_x) * u_y + v_xyAlongu);
        }
      }
      else if constexpr (DerivativeOrder
          == 3) {  // dd(dot(u,v))/(dxdydz) =  u_{x,y,z} * v + u_{x,y} * v_z + u_{x,z}*v_y + u_x*v_{y,z} + u_{y,z}* v_x + u_y* v_{x,z} + u_z * v_{x,y} + u * v_{x,y,z}
        //localFunctionArgs is here argsForDxyz
        const auto u = evaluateFunctionImpl(this->l, localFunctionArgs).getValue();
        const auto v = evaluateFunctionImpl(this->r, localFunctionArgs).getValue();

        const auto argsForDx  = localFunctionArgs.extractSpatialOrFirstWrtArg();
        const auto argsForDyz = localFunctionArgs.extractSecondWrtArgOrFirstNonSpatial();

        const auto [argsForDy, argsForDz] = localFunctionArgs.extractWrtArgsTwoCoeffsToSingleCoeff();

        const auto& [u_x,u_y,u_z] = evaluateFirstOrderDerivativesImpl(this->l, localFunctionArgs);
        const auto& [v_x,v_y,v_z] = evaluateFirstOrderDerivativesImpl(this->r, localFunctionArgs);
        const auto& [u_xy,u_xz] = evaluateSecondOrderDerivativesImpl(this->l, localFunctionArgs);

//        const auto argsForDxy = joinWRTArgs(argsForDx,argsForDy);
//        const auto argsForDxz = joinWRTArgs(argsForDx,argsForDz);

        const auto alonguArgs = createWithAlong(localFunctionArgs, along(v));
        const auto alongvArgs = createWithAlong(localFunctionArgs, along(u));
        const auto argsForDyzalongv_xArgs = createWithAlong(argsForDyz, along(v_x));
        const auto argsForDyzalongu_xArgs = createWithAlong(argsForDyz, along(u_x));

        const auto u_xyzAlongv = evaluateDerivativeImpl(this->l, alongvArgs);
        const auto v_xyzAlongu = evaluateDerivativeImpl(this->r, alonguArgs);
        const auto u_yzAlongvx = evaluateDerivativeImpl(this->l, argsForDyzalongv_xArgs);
        const auto v_yzAlongux        = evaluateDerivativeImpl(this->r, argsForDyzalongu_xArgs);
//        const auto u_x        = evaluateDerivativeImpl(this->l, argsForDx);
//        const auto u_y        = evaluateDerivativeImpl(this->l, argsForDy);
//        const auto u_z        = evaluateDerivativeImpl(this->l, argsForDz);

//        const auto v_x        = evaluateDerivativeImpl(this->r, argsForDx);
//        const auto v_y        = evaluateDerivativeImpl(this->r, argsForDy);
//        const auto v_z        = evaluateDerivativeImpl(this->r, argsForDz);


//        const auto u_xy = evaluateDerivativeImpl(this->l, argsForDxy);
//        const auto u_xz = evaluateDerivativeImpl(this->l, argsForDxz);


//        const auto v_xy        = evaluateDerivativeImpl(this->r, argsForDxy);
//        const auto v_xz        = evaluateDerivativeImpl(this->r, argsForDxz);

        return eval(u_xyzAlongv + u_xy*v_z + u_xz*v_y + v_yzAlongux + u_yzAlongvx + v_xz*u_y  + v_xy*u_z  + v_xyzAlongu);
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