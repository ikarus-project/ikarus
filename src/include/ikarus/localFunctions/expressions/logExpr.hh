//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
namespace Ikarus {

  template <typename E1>
  class LogExpr : public UnaryLocalFunctionExpression<LogExpr, E1> {
  public:
    using Base = UnaryLocalFunctionExpression<LogExpr, E1>;
    using Base::UnaryLocalFunctionExpression;
    using Traits                   = LocalFunctionTraits<LogExpr>;
    static constexpr int valueSize = 1;
    static constexpr int gridDim   = Traits::gridDim;
    using ctype                    = typename Traits::ctype;

    using E1Raw = std::remove_cvref_t<E1>;

    template <size_t ID_ = 0>
    static constexpr int orderID = nonLinear;

    template <typename LFArgs>
    auto evaluateValueOfExpression(const LFArgs& lfArgs) const {
      return Eigen::Vector<ctype, 1>(Ikarus::eval(log(evaluateFunctionImpl(this->m(), lfArgs)[0])));
    }

    template <int DerivativeOrder, typename LFArgs>
    auto evaluateDerivativeOfExpression(const LFArgs& lfArgs) const {
      const auto u = evaluateFunctionImpl(this->m(), lfArgs)[0];
      if constexpr (DerivativeOrder == 1)
      {
        const auto u_x = evaluateDerivativeImpl(this->m(), lfArgs);
        return Ikarus::eval(u_x / u);
      }
      else if constexpr (DerivativeOrder == 2)
      {
        const auto& [u_x, u_y]    = evaluateFirstOrderDerivativesImpl(this->m(), lfArgs);
        const auto u_xy           = evaluateDerivativeImpl(this->m(), lfArgs);

        if constexpr (LFArgs::hasOneSpatialAll and LFArgs::hasSingleCoeff) {
          std::array<std::remove_cvref_t<decltype(Ikarus::eval(u_x.col(0) * u_y))>, gridDim> res;
          for (int i = 0; i < gridDim; ++i)
            res[i] = Ikarus::eval(-u_x.col(i) * u_y /(u*u) + u_xy[i]/u);
          return res;
        } else {  // one single spatial and one coeff derivative OR two coeffs
          return Ikarus::eval(-transpose(u_x) * u_y /(u*u) + u_xy/u);
        }
      } else if constexpr (DerivativeOrder
                           == 3) {
        const auto& [u_x, u_y, u_z] = evaluateFirstOrderDerivativesImpl(this->m(), lfArgs);
        const auto& [u_xy, u_xz]    = evaluateSecondOrderDerivativesImpl(this->m(), lfArgs);

        const auto u_xyz = evaluateDerivativeImpl(this->m(), lfArgs);

        const auto argsForDyz = lfArgs.extractSecondWrtArgOrFirstNonSpatial();
        const auto u_yz       = evaluateDerivativeImpl(this->m(), argsForDyz);

        if constexpr (LFArgs::hasOneSpatialSingle) {
          static_assert(decltype(u_x)::ColsAtCompileTime == 1);
          static_assert(decltype(u_x)::RowsAtCompileTime == 1);

          return eval(u_xyz/u
                      - (transpose(u_xy)*u_z + transpose(u_y)*u_xz + transpose(u_yz)*u_x[0])/(u*u)
                      + 2*u_x[0]* transpose(u_y)*u_z/(u*u*u));

        } else if constexpr (LFArgs::hasOneSpatialAll) {
          const auto& alongMatrix = std::get<0>(lfArgs.alongArgs.args);
          std::remove_cvref_t<decltype(Ikarus::eval(transpose(u_y) * u_z))> res;
          res = u_xyz / u;

          for (int i = 0; i < gridDim; ++i)
            res += alongMatrix(0, i) *
                   (- (transpose(u_xy[i])*u_z + transpose(u_y)*u_xz[i] + transpose(u_yz)*u_x[i])/(u*u)
                    + 2*u_x[i]* transpose(u_y)*u_z/(u*u*u));
//            res += alongMatrix(0, i)
//                   * (3 * u_x[i] * transpose(u_y) * u_z / (8 * pow(u, 5.0 / 2.0))
//                      - (transpose(u_y) * u_xz[i] + u_x[i] * transpose(u_yz) + transpose(u_xy[i]) * u_z)
//                            * udiv4timesppow32);
          return res;
        }
      }
    }
  };

  template <typename E1>
  struct LocalFunctionTraits<LogExpr<E1>> : public LocalFunctionTraits<std::remove_cvref_t<E1>> {};

  template <typename E1>
  requires IsLocalFunction<E1>
  constexpr auto log(E1&& u) {
    static_assert(std::remove_cvref_t<E1>::valueSize == 1,
                  "Log expression only defined for scalar valued local functions.");
    return LogExpr<E1>(std::forward<E1>(u));
  }

}  // namespace Ikarus