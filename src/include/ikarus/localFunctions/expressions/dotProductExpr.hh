//
// Created by lex on 4/25/22.
//

#pragma once
#include "rebind.hh"

#include <ikarus/localFunctions/expressions/binaryExpr.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
namespace Ikarus {

  template <typename E1, typename E2>
  class LocalFunctionDot : public BinaryLocalFunctionExpression<LocalFunctionDot, E1, E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionDot, E1, E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionDot>;
    /** \brief Type used for coordinates */
    using ctype                    = typename Traits::ctype;
    static constexpr int valueSize = Traits::valueSize;
    static constexpr int gridDim   = Traits::gridDim;

    template <size_t ID_ = 0>
    static constexpr int orderID
        = std::min(Base::E1Raw::template order<ID_>() + Base::E2Raw::template order<ID_>(), nonLinear);

    template <typename LFArgs>
    auto evaluateValueOfExpression(const LFArgs &lfArgs) const {
      const auto u = evaluateFunctionImpl(this->l(), lfArgs);
      const auto v = evaluateFunctionImpl(this->r(), lfArgs);
      return Eigen::Vector<ctype, 1>((u.transpose() * v).trace());
    }

    template <int DerivativeOrder, typename LFArgs>
    auto evaluateDerivativeOfExpression(const LFArgs &lfArgs) const {
      const auto u = evaluateFunctionImpl(this->l(), lfArgs);
      const auto v = evaluateFunctionImpl(this->r(), lfArgs);
      if constexpr (DerivativeOrder == 1)  // d(dot(u,v))/dx =  u_x * v+ u*v_x
      {
        const auto u_x = evaluateDerivativeImpl(this->l(), lfArgs);
        const auto v_x = evaluateDerivativeImpl(this->r(), lfArgs);
        return Ikarus::eval(v.transpose() * u_x + u.transpose() * v_x);
      } else if constexpr (DerivativeOrder
                           == 2) {  // dd(dot(u,v))/(dxdy) =  u_{x,y} * v + u_x*v_y + u_y* v_x + u * v_{x,y}
        const auto &[u_x, u_y] = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs);
        const auto &[v_x, v_y] = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs);
        if constexpr (LFArgs::hasNoSpatial and LFArgs::hasTwoCoeff) {
          const auto alonguArgs = replaceAlong(lfArgs, along(v));
          const auto alongvArgs = replaceAlong(lfArgs, along(u));

          const auto u_xyAlongv = evaluateDerivativeImpl(this->l(), alongvArgs);
          const auto v_xyAlongu = evaluateDerivativeImpl(this->r(), alonguArgs);

          return Ikarus::eval(u_xyAlongv + transpose(u_x) * v_y + transpose(v_x) * u_y + v_xyAlongu);
        } else if constexpr (LFArgs::hasOneSpatial and LFArgs::hasSingleCoeff) {
          const auto u_xy = evaluateDerivativeImpl(this->l(), lfArgs);
          const auto v_xy = evaluateDerivativeImpl(this->r(), lfArgs);
          if constexpr (LFArgs::hasOneSpatialSingle and LFArgs::hasSingleCoeff) {
            return Ikarus::eval(transpose(v) * u_xy + transpose(u_x) * v_y + transpose(v_x) * u_y
                                + transpose(u) * v_xy);
          } else if constexpr (LFArgs::hasOneSpatialAll and LFArgs::hasSingleCoeff) {
            std::array<std::remove_cvref_t<decltype(Ikarus::eval(transpose(v) * u_xy[0]))>, gridDim> res;
            for (int i = 0; i < gridDim; ++i)
              res[i] = Ikarus::eval(transpose(v) * u_xy[i] + transpose(u_x.col(i)) * v_y + transpose(v_x.col(i)) * u_y
                                    + transpose(u) * v_xy[i]);
            return res;
          }
        }
      } else if constexpr (DerivativeOrder
                           == 3) {  // dd(dot(u,v))/(dxdydz) =  u_{x,y,z} * v + u_{x,y} * v_z + u_{x,z}*v_y +
                                    // u_x*v_{y,z} + u_{y,z}* v_x + u_y* v_{x,z} + u_z * v_{x,y} + u * v_{x,y,z}
        if constexpr (LFArgs::hasOneSpatialSingle) {
          const auto argsForDyz = lfArgs.extractSecondWrtArgOrFirstNonSpatial();

          const auto &[u_x, u_y, u_z] = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs);
          const auto &[v_x, v_y, v_z] = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs);
          const auto &[u_xy, u_xz]    = evaluateSecondOrderDerivativesImpl(this->l(), lfArgs);
          const auto &[v_xy, v_xz]    = evaluateSecondOrderDerivativesImpl(this->r(), lfArgs);

          const auto alonguArgs             = replaceAlong(lfArgs, along(u));
          const auto alongvArgs             = replaceAlong(lfArgs, along(v));
          const auto argsForDyzalongv_xArgs = replaceAlong(argsForDyz, along(v_x));
          const auto argsForDyzalongu_xArgs = replaceAlong(argsForDyz, along(u_x));

          const auto u_xyzAlongv = evaluateDerivativeImpl(this->l(), alongvArgs);
          const auto v_xyzAlongu = evaluateDerivativeImpl(this->r(), alonguArgs);
          const auto u_yzAlongvx = evaluateDerivativeImpl(this->l(), argsForDyzalongv_xArgs);
          const auto v_yzAlongux = evaluateDerivativeImpl(this->r(), argsForDyzalongu_xArgs);

          return Ikarus::eval(u_xyzAlongv + transpose(u_xy) * v_z + transpose(u_xz) * v_y + v_yzAlongux + u_yzAlongvx
                              + transpose(v_xz) * u_y + transpose(v_xy) * u_z + v_xyzAlongu);
        } else if constexpr (LFArgs::hasOneSpatialAll) {
          // check that the along argument has the correct size
          const auto &alongMatrix = std::get<0>(lfArgs.alongArgs.args);
          static_assert(alongMatrix.ColsAtCompileTime == gridDim);
          static_assert(alongMatrix.RowsAtCompileTime == 1);

          const auto uTimesA = eval(u * alongMatrix);
          const auto vTimesA = eval(v * alongMatrix);
          static_assert(uTimesA.RowsAtCompileTime == Base::E1Raw::valueSize);
          static_assert(vTimesA.RowsAtCompileTime == Base::E2Raw::valueSize);
          static_assert(uTimesA.ColsAtCompileTime == gridDim);
          static_assert(vTimesA.ColsAtCompileTime == gridDim);

          const auto &[gradu, u_c0, u_c1]  = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs);
          const auto &[gradv, v_c0, v_c1]  = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs);
          const auto &[gradu_c0, gradu_c1] = evaluateSecondOrderDerivativesImpl(this->l(), lfArgs);
          const auto &[gradv_c0, gradv_c1] = evaluateSecondOrderDerivativesImpl(this->r(), lfArgs);

          const auto graduTimesA = (gradu * alongMatrix.transpose()).eval();
          const auto gradvTimesA = (gradv * alongMatrix.transpose()).eval();
          static_assert(graduTimesA.RowsAtCompileTime == Base::E1Raw::valueSize);
          static_assert(graduTimesA.ColsAtCompileTime == 1);
          static_assert(gradvTimesA.RowsAtCompileTime == Base::E2Raw::valueSize);
          static_assert(gradvTimesA.ColsAtCompileTime == 1);

          const auto argsForDyz = lfArgs.extractSecondWrtArgOrFirstNonSpatial();

          const auto alonguAArgs          = replaceAlong(lfArgs, along(uTimesA));
          const auto alongvAArgs          = replaceAlong(lfArgs, along(vTimesA));
          const auto alonggraduTimesAArgs = replaceAlong(argsForDyz, along(graduTimesA));
          const auto alonggradvTimesAArgs = replaceAlong(argsForDyz, along(gradvTimesA));

          const auto u_xyzAlongv            = evaluateDerivativeImpl(this->l(), alongvAArgs);
          const auto v_xyzAlongu            = evaluateDerivativeImpl(this->r(), alonguAArgs);
          const auto v_c0c1AlongGraduTimesA = evaluateDerivativeImpl(this->r(), alonggraduTimesAArgs);
          const auto u_c0c1AlongGradvTimesA = evaluateDerivativeImpl(this->l(), alonggradvTimesAArgs);
          decltype(eval(u_xyzAlongv)) res;

          res = u_xyzAlongv + v_xyzAlongu + v_c0c1AlongGraduTimesA + u_c0c1AlongGradvTimesA;
          for (int i = 0; i < gridDim; ++i)
            res += (transpose(u_c1) * gradv_c0[i] + transpose(v_c1) * gradu_c0[i] + transpose(v_c0) * gradu_c1[i]
                    + transpose(u_c0) * gradv_c1[i])
                   * alongMatrix(0, i);

          return res;

        } else
          static_assert(
              LFArgs::hasOneSpatialSingle or LFArgs::hasOneSpatialAll,
              "Only a spatial single direction or all spatial directions are supported. You should not end up here.");
      } else
        static_assert(DerivativeOrder > 3 or DerivativeOrder < 1,
                      "Only first, second and third order derivatives are supported.");
    }
  };

  template <typename E1, typename E2>
  struct LocalFunctionTraits<LocalFunctionDot<E1, E2>> {
    using E1Raw = std::remove_cvref_t<E1>;
    using E2Raw = std::remove_cvref_t<E2>;
    /** \brief Size of the function value */
    static constexpr int valueSize = 1;
    /** \brief Type for the points for evaluation, usually the integration points */
    using DomainType = std::common_type_t<typename E1Raw::DomainType, typename E2Raw::DomainType>;
    /** \brief Type used for coordinates */
    using ctype = std::common_type_t<typename E1Raw::ctype, typename E1Raw::ctype>;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = E1Raw::gridDim;
  };

  template <typename E1, typename E2>
  requires IsLocalFunction<E1, E2>
  constexpr auto dot(E1 &&u, E2 &&v) { return LocalFunctionDot<E1, E2>(std::forward<E1>(u), std::forward<E2>(v)); }

}  // namespace Ikarus