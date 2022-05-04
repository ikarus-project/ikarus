//
// Created by lex on 4/25/22.
//

#pragma once
#include <ikarus/localFunctions/expressions/binaryExpr.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/localFunctions/expressions/rebind.hh>

namespace Ikarus {

  template <typename E1, typename E2>
  class LocalFunctionSum : public BinaryLocalFunctionExpression<LocalFunctionSum, E1, E2> {
  public:
    using Base = BinaryLocalFunctionExpression<LocalFunctionSum, E1, E2>;
    using Base::BinaryLocalFunctionExpression;
    using Traits = LocalFunctionTraits<LocalFunctionSum>;

    template<size_t ID_=0>
    static constexpr int order = std::max(std::remove_cvref_t<E1>::template order<ID_>,std::remove_cvref_t<E2>::template order<ID_>)  ;

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
      return Ikarus::eval(evaluateFunctionImpl(this->l(), localFunctionArgs)
                  + evaluateFunctionImpl(this->r(), localFunctionArgs));
    }

    template <int DerivativeOrder, typename LocalFunctionEvaluationArgs_>
    auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return Ikarus::eval(evaluateDerivativeImpl(this->l(), localFunctionArgs)
                  + evaluateDerivativeImpl(this->r(), localFunctionArgs));
    }
  };

  template <typename E1, typename E2>
  struct LocalFunctionTraits<LocalFunctionSum<E1, E2>> : public LocalFunctionTraits<std::remove_cvref_t<E1>> {
    using Base = LocalFunctionTraits<std::remove_cvref_t<E1>>;
  };

  template <typename E1, typename E2> requires IsLocalFunction<E1,E2>
  constexpr auto operator+(E1 && u, E2 && v) {
    static_assert(Concepts::AddAble<typename std::remove_cvref_t<E1>::FunctionReturnType, typename std::remove_cvref_t<E2>::FunctionReturnType>,
                  "The function values of your local functions are not addable!");
    return LocalFunctionSum<E1,E2>(std::forward<E1>(u), std::forward<E2>(v));
  }
}  // namespace Ikarus