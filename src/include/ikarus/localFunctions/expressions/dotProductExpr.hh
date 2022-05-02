//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/binaryExpr.hh>
#include <ikarus/manifolds/realTuple.hh>
#include "rebind.hh"
namespace Ikarus {

  template <typename E1, typename E2>
  class LocalFunctionDot : public BinaryLocalFunctionExpression<LocalFunctionDot, E1, E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionDot, E1, E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionDot>;
    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    static constexpr int valueSize =  Traits::valueSize;

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      const auto u = evaluateFunctionImpl(this->l(), localFunctionArgs).getValue();
      const auto v = evaluateFunctionImpl(this->r(), localFunctionArgs).getValue();
      return Ikarus::RealTuple<ctype,1>( Eigen::Matrix<ctype,1,1>(u.dot(v)));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      const auto u = evaluateFunctionImpl(this->l(), localFunctionArgs).getValue();
      const auto v = evaluateFunctionImpl(this->r(), localFunctionArgs).getValue();
      if constexpr (DerivativeOrder == 1)  // d(dot(u,v))/dx =  u_x * v+ u*v_x
      {
        const auto u_x = evaluateDerivativeImpl(this->l(), localFunctionArgs);
        const auto v_x = evaluateDerivativeImpl(this->r(), localFunctionArgs);
        return Ikarus::eval(v.transpose()*u_x + u.transpose()*v_x);
      }
      else if constexpr (DerivativeOrder == 2 ) {  // dd(dot(u,v))/(dxdy) =  u_{x,y} * v + u_x*v_y + u_y* v_x + u * v_{x,y}
        const auto& [u_x,u_y] = evaluateFirstOrderDerivativesImpl(this->l(), localFunctionArgs);
        const auto& [v_x,v_y] = evaluateFirstOrderDerivativesImpl(this->r(), localFunctionArgs);
        if constexpr (LocalFunctionEvaluationArgs_::hasOneSpatialSingle
                      and LocalFunctionEvaluationArgs_::hasSingleCoeff) {

          const auto u_xy = evaluateDerivativeImpl(this->l(), localFunctionArgs);
          const auto v_xy = evaluateDerivativeImpl(this->r(), localFunctionArgs);

          return Ikarus::eval(transpose(v) * u_xy + transpose(u_x) * v_y + transpose(v_x) * u_y + transpose(u) * v_xy);
        } else if constexpr (LocalFunctionEvaluationArgs_::hasNoSpatial and LocalFunctionEvaluationArgs_::hasTwoCoeff) {
          const auto alonguArgs = addAlong(localFunctionArgs, along(v));
          const auto alongvArgs = addAlong(localFunctionArgs, along(u));

          const auto u_xyAlongv = evaluateDerivativeImpl(this->l(), alongvArgs);
          const auto v_xyAlongu = evaluateDerivativeImpl(this->r(), alonguArgs);

          return Ikarus::eval(u_xyAlongv + transpose(u_x) * v_y + transpose(v_x) * u_y + v_xyAlongu);
        }
      }
      else if constexpr (DerivativeOrder
          == 3) {  // dd(dot(u,v))/(dxdydz) =  u_{x,y,z} * v + u_{x,y} * v_z + u_{x,z}*v_y + u_x*v_{y,z} + u_{y,z}* v_x + u_y* v_{x,z} + u_z * v_{x,y} + u * v_{x,y,z}
        //localFunctionArgs is here argsForDxyz

        const auto argsForDyz = localFunctionArgs.extractSecondWrtArgOrFirstNonSpatial();

        const auto& [u_x,u_y,u_z] = evaluateFirstOrderDerivativesImpl(this->l(), localFunctionArgs);
        const auto& [v_x,v_y,v_z] = evaluateFirstOrderDerivativesImpl(this->r(), localFunctionArgs);
        const auto& [u_xy,u_xz] = evaluateSecondOrderDerivativesImpl(this->l(), localFunctionArgs);
        const auto& [v_xy,v_xz] = evaluateSecondOrderDerivativesImpl(this->r(), localFunctionArgs);

        const auto alonguArgs = addAlong(localFunctionArgs, along(v));
        const auto alongvArgs = addAlong(localFunctionArgs, along(u));
        const auto argsForDyzalongv_xArgs = addAlong(argsForDyz, along(v_x));
        const auto argsForDyzalongu_xArgs = addAlong(argsForDyz, along(u_x));

        const auto u_xyzAlongv = evaluateDerivativeImpl(this->l(), alongvArgs);
        const auto v_xyzAlongu = evaluateDerivativeImpl(this->r(), alonguArgs);
        const auto u_yzAlongvx = evaluateDerivativeImpl(this->l(), argsForDyzalongv_xArgs);
        const auto v_yzAlongux = evaluateDerivativeImpl(this->r(), argsForDyzalongu_xArgs);

        return Ikarus::eval(u_xyzAlongv + transpose(u_xy)*v_z + transpose(u_xz)*v_y + v_yzAlongux + u_yzAlongvx + transpose(v_xz)*u_y  + transpose(v_xy)*u_z  + v_xyzAlongu);
      }else
        static_assert(DerivativeOrder > 3 or DerivativeOrder<1,
                      "Only first, secondand third order derivatives are supported.");
    }

//    auto clone()const{
//      return LocalFunctionDot(this->l().clone(), this->r().clone());
//    }
  };

  template <typename E1, typename E2>
  struct LocalFunctionTraits<LocalFunctionDot<E1, E2>> {
    using E1Raw = std::remove_cvref_t<E1>;
    /** \brief Size of the function value */
    static constexpr int valueSize = 1;
    /** \brief Type for the points for evaluation, usually the integration points */
    using DomainType = typename E1Raw::DomainType;
    /** \brief Type used for coordinates */
    using ctype = typename E1Raw::ctype;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = E1Raw::gridDim;
  };

  template <typename E1, typename E2> requires IsLocalFunction<E1,E2>
  constexpr auto dot(E1 && u, E2 && v) {
    return LocalFunctionDot<E1,E2>(std::forward<E1>(u), std::forward<E2>(v));
  }

}  // namespace Ikarus