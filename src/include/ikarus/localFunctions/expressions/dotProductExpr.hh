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

    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    //    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = Traits::valueSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Traits::gridDim;
    /** \brief Type for coordinate vector in world space */
    using FunctionReturnType = typename Traits::FunctionReturnType;

    template <typename LocalFunctionEvaluationArgs_>
    auto evaluateValueOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return evaluateFunctionImpl(this->_u, localFunctionArgs)
          .getValue()
          .dot(evaluateFunctionImpl(this->_v, localFunctionArgs).getValue());
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      static_assert(DerivativeOrder == 1, "For testing only first order derivatives are supported.");
      if constexpr (DerivativeOrder == 1)
        return eval((transpose(evaluateDerivativeImpl(this->_u, localFunctionArgs))
                      * evaluateFunctionImpl(this->_v,localFunctionArgs).getValue()).transpose()
                  + evaluateFunctionImpl(this->_u,localFunctionArgs).getValue().transpose()
                        * evaluateDerivativeImpl(this->_v, localFunctionArgs));
      else if constexpr (DerivativeOrder == 2) {
        const auto args = LocalFunctionEvaluationArgs(localFunctionArgs,along(evaluateFunctionImpl(this->_v, localFunctionArgs).getValue()));
        const auto argsFirstOrderDeriv = LocalFunctionEvaluationArgs(localFunctionArgs,wrt(std::));
//
//        return eval((evaluateDerivativeImpl(this->_u, args).transpose()
//                     * evaluateFunctionImpl(this->_v, localFunctionArgs).getValue())
//                        .transpose()
//                    + evaluateFunctionImpl(this->_u, localFunctionArgs).getValue().transpose()
//                          * evaluateDerivativeImpl(this->_v, localFunctionArgs));
      }
    }
  };

  template <typename E1, typename E2>
  struct LocalFunctionTraits<LocalFunctionDot<E1, E2>> {
    /** \brief Type used for coordinates */
    using ctype                    = typename E1::ctype;
    static constexpr int valueSize = 1;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = E1::gridDim;
    /** \brief Type for the return value */
    using FunctionReturnType = ctype;
    /** \brief Type for the derivatives wrt. the coeffiecients */

    /** \brief Type for the points for evaluation, usually the integration points */
    using DomainType = typename E1::DomainType;
  };

  template <typename E1, typename E2>
  LocalFunctionDot<E1, E2> dot(LocalFunctionInterface<E1> const& u, LocalFunctionInterface<E2> const& v) {
    return LocalFunctionDot<E1, E2>(u, v);
  }

}  // namespace Ikarus