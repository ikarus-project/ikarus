//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
namespace Ikarus {

  template <typename E1>
  class SqrtExpr : public UnaryLocalFunctionExpression<SqrtExpr, E1> {
  public:
    using Base = UnaryLocalFunctionExpression<SqrtExpr, E1>;
    using Base::UnaryLocalFunctionExpression;
    using Traits                   = LocalFunctionTraits<SqrtExpr>;
    static constexpr int valueSize = 1;
    static constexpr int gridDim   = Traits::gridDim;
    using ctype                    = typename Traits::ctype;

    using E1Raw = std::remove_cvref_t<E1>;

    template <size_t ID_ = 0>
    static constexpr int orderID = nonLinear;

    template <typename LFArgs>
    auto evaluateValueOfExpression(const LFArgs& lfArgs) const {
      return Eigen::Vector<ctype, 1>(Ikarus::eval(sqrt(evaluateFunctionImpl(this->m(), lfArgs)[0])));
    }

    template <int DerivativeOrder, typename LFArgs>
    auto evaluateDerivativeOfExpression(const LFArgs& lfArgs) const {
      const auto u = evaluateFunctionImpl(this->m(), lfArgs)[0];
      if constexpr (DerivativeOrder == 1)  // d(sqrt(u(x)))/(dx) =  u_x /(2*sqrt(u(x))
      {
        const auto u_x = evaluateDerivativeImpl(this->m(), lfArgs);
        return Ikarus::eval(u_x / (2 * sqrt(u)));
      } else if constexpr (DerivativeOrder == 2) {  // d^3(qrt(u(x,y)))/(dxdy) = - u_x*u_y/(4*u^(3/2)) +
                                                    // u_xy/(2*sqrt(u))
        const auto& [u_x, u_y]    = evaluateFirstOrderDerivativesImpl(this->m(), lfArgs);
        const auto u_xy           = evaluateDerivativeImpl(this->m(), lfArgs);
        const auto u_yTimesfactor = Ikarus::eval(u_y / (4 * std::pow(u, 3.0 / 2.0)));
        if constexpr (LFArgs::hasOneSpatialAll and LFArgs::hasSingleCoeff) {
          std::array<std::remove_cvref_t<decltype(Ikarus::eval(u_x.col(0) * u_y))>, gridDim> res;
          for (int i = 0; i < gridDim; ++i)
            res[i] = Ikarus::eval(-u_x.col(i) * u_yTimesfactor + u_xy[i] * 1.0 / (2 * sqrt(u)));
          return res;
        } else {  // one spatial and one coeff derivative
          return Ikarus::eval(-transpose(u_x) * u_yTimesfactor + u_xy / (2 * sqrt(u)));
        }
      } else if constexpr (DerivativeOrder
                           == 3) {  // d^3(qrt(u(x,y,z)))/(dxdydz) =(3*u_x*u_y*u_z)/(8*u^(5/2)) - u_xz*u_y/(4*u^(3/2)) -
                                    // u_x*u_yz/(4*u^(3/2)) - u_xy*u_z/(4*u^(3/2)) + u_xyz/(2*sqrt(u))
        const auto& [u_x, u_y, u_z] = evaluateFirstOrderDerivativesImpl(this->m(), lfArgs);
        const auto& [u_xy, u_xz]    = evaluateSecondOrderDerivativesImpl(this->m(), lfArgs);

        const auto u_xyz = evaluateDerivativeImpl(this->m(), lfArgs);

        const auto argsForDyz = lfArgs.extractSecondWrtArgOrFirstNonSpatial();
        const auto u_yz       = evaluateDerivativeImpl(this->m(), argsForDyz);

        const ctype udiv4timesppow32 = 1.0 / (4.0 * pow(u, 3.0 / 2.0));
        if constexpr (LFArgs::hasOneSpatialSingle) {
          static_assert(decltype(u_x)::ColsAtCompileTime == 1);
          static_assert(decltype(u_x)::RowsAtCompileTime == 1);

          return eval(u_xyz / (2.0 * sqrt(u)) + (3.0 * u_x[0] * transpose(u_y) * u_z) / (8.0 * pow(u, 5.0 / 2.0))
                      - ((transpose(u_y) * u_xz + u_x[0] * transpose(u_yz) + transpose(u_xy) * u_z))
                            * udiv4timesppow32);
        } else if constexpr (LFArgs::hasOneSpatialAll) {
          const auto& alongMatrix = std::get<0>(lfArgs.alongArgs.args);
          std::remove_cvref_t<decltype(Ikarus::eval(transpose(u_y) * u_z))> res;
          res = u_xyz / (2.0 * sqrt(u));

          for (int i = 0; i < gridDim; ++i)
            res += alongMatrix(0, i)
                   * (3 * u_x[i] * transpose(u_y) * u_z / (8 * pow(u, 5.0 / 2.0))
                      - (transpose(u_y) * u_xz[i] + u_x[i] * transpose(u_yz) + transpose(u_xy[i]) * u_z)
                            * udiv4timesppow32);
          return res;
        }
      }
    }
  };

  template <typename E1>
  struct LocalFunctionTraits<SqrtExpr<E1>> : public LocalFunctionTraits<std::remove_cvref_t<E1>> {};

  template <typename E1>
    requires IsLocalFunction<E1>
  constexpr auto sqrt(E1&& u) {
    static_assert(std::remove_cvref_t<E1>::valueSize == 1,
                  "Sqrt expression only defined for scalar valued local functions.");
    return SqrtExpr<E1>(std::forward<E1>(u));
  }

}  // namespace Ikarus