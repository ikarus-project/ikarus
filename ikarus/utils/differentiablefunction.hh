// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file differentiablefunction.hh
 * \brief Provides a DifferentiableFunction class for handling differentiable Functions.
 *
 */

#pragma once
#include <tuple>
#include <type_traits>

#include <dune/common/hybridutilities.hh>
#include <dune/functions/common/differentiablefunctionfromcallables.hh>

#include <ikarus/utils/derivativetraits.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

namespace Impl {

  /**
   * \brief Represents a tuple of functions.
   *
   * \tparam Args The argument types.
   */
  template <typename... Args>
  struct Functions
  {
    std::tuple<Args...> args;
  };
} // namespace Impl

/**
 * \brief Creates a Functions object.
 *
 * \tparam Args The argument types.
 * \param args The tuple of arguments.
 * \return auto The Functions object.
 */
template <typename... Args>
auto functions(Args&&... args) {
  return Impl::Functions<std::remove_cvref_t<Args>...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

#ifndef DOXYGEN
template <class Signature, template <class> class DerivativeTraits, class... F>
class DifferentiableFunction;
#endif

/**
 * \brief DifferentiableFunction is a class taking several callables.
 * The function are assumed to be derivatives of each other w.r.t. the argument
 *
 * \tparam DerivativeArgs The types of derivative arguments.
 * \tparam ParameterArgs The types of parameter arguments.
 */
template <class Range, class D, template <class> class DerivativeTraits, class F, class... FF>
class DifferentiableFunction<Range(D), DerivativeTraits, F, FF...>
    : private Dune::Functions::DifferentiableFunctionFromCallables<Range(D), DerivativeTraits, F, FF...>
{
  using Base = Dune::Functions::DifferentiableFunctionFromCallables<Range(D), DerivativeTraits, F, FF...>;

public:
  using Domain                       = std::remove_cvref_t<D>;
  static constexpr auto nDerivatives = sizeof...(FF);

  template <class... FFF>
  DifferentiableFunction(FFF&&... f)
      : Base(std::forward<FFF>(f)...) {}
  using Traits = DerivativeTraitsFromCallables<Impl::Functions<F, FF...>, D>;

  using Derivative = DifferentiableFunction<typename Traits::template Range<1>(D), DerivativeTraits, FF...>;

  Range operator()(const D& x) const { return Base::operator()(x); }

  /**
   * \brief Get derivative of DifferentiableFunction
   *
   */
  friend Derivative derivative(const DifferentiableFunction& t) {
    auto df = derivative(static_cast<const Base&>(t));
    return Derivative(df);
  }
};

/**
 * \brief Factory method for DifferentiableFunction
 * It is a function taking several callables and the argument of these functions to derive the correct signatures for
 * the functions (which could be templated lambdas) The function are assumed to be derivatives of each other w.r.t. the
 * parameter
 *
 * \tparam F The types of  functions.
 * \tparam Arg The types of the argument.
 */
template <typename... F, typename Arg>
auto makeDifferentiableFunction(const Impl::Functions<F...>& derivativesFunctions, const Arg& parameter) {
  DerivativeTraitsFromCallables t(derivativesFunctions, parameter);
  using DerivTraits = decltype(t);
  auto la           = []<typename... FF>(FF&&... f) {
    return Ikarus::DifferentiableFunction<typename DerivTraits::template Range<0>(typename DerivTraits::Domain),
                                                    DerivTraits::template DerivativeTraits, std::remove_cvref_t<FF>...>(
        std::forward<FF>(f)...);
  };
  return std::apply(la, derivativesFunctions.args);
}

} // namespace Ikarus
