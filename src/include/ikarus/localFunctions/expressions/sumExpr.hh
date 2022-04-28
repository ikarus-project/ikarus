//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/binaryExpr.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename E1, typename E2>
  class LocalFunctionSum : public BinaryLocalFunctionExpression<LocalFunctionSum<E1, E2>, E1, E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionSum<E1, E2>, E1, E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionSum>;

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
      return eval(evaluateFunctionImpl(this->l, localFunctionArgs)
                  + evaluateFunctionImpl(this->r, localFunctionArgs));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return eval(evaluateDerivativeImpl(this->l, localFunctionArgs)
                  + evaluateDerivativeImpl(this->r, localFunctionArgs));
    }
  };

  template <typename E1, typename E2>
  struct LocalFunctionTraits<LocalFunctionSum<E1, E2>> : public LocalFunctionTraits<E1> {
    using Base = LocalFunctionTraits<E1>;
  };

  template <typename E1, typename E2>
  constexpr LocalFunctionSum<E1, E2> operator+(LocalFunctionInterface<E1> const& u, LocalFunctionInterface<E2> const& v) {
    static_assert(Concepts::AddAble<typename E1::FunctionReturnType, typename E2::FunctionReturnType>,
                  "The function values of your local functions are not addable!");

    return LocalFunctionSum<E1, E2>(u, v);
  }
}  // namespace Ikarus