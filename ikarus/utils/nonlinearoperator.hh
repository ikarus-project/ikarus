// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearoperator.hh
 * \brief Provides a NonLinearOperator class for handling nonlinear operators.
 *
 */

#pragma once
#include <functional>

#include <dune/common/hybridutilities.hh>

#include <ikarus/utils/traits.hh>

namespace Ikarus {

namespace Impl {
  /**
   * \brief Helper function to apply a function and remove reference wrappers.
   *
   * \tparam F The function type.
   * \tparam Tuple The tuple type.
   * \tparam I Index sequence for tuple elements.
   * \param f The function to be applied.
   * \param t The tuple of arguments.
   * \return constexpr decltype(auto) The result after applying the function and removing reference wrappers.
   */
  template <class F, class Tuple, std::size_t... I>
  constexpr decltype(auto) applyAndRemoveRefererenceWrapper(F&& f, Tuple&& t, std::index_sequence<I...>) {
    return std::invoke(std::forward<F>(f),
                       std::get<I>(std::forward<Tuple>(t)).get()...); //.get gets the impl type of std::referenceWrapper
  }

  /**
   * \brief Helper function to apply a function and remove reference wrappers.
   *
   * \tparam F The function type.
   * \tparam Tuple The tuple type.
   * \param f The function to be applied.
   * \param t The tuple of arguments.
   * \return constexpr decltype(auto) The result after applying the function and removing reference wrappers.
   */
  template <class F, class Tuple>
  constexpr decltype(auto) applyAndRemoveReferenceWrapper(F&& f, Tuple&& t) {
    return applyAndRemoveRefererenceWrapper(
        std::forward<F>(f), std::forward<Tuple>(t),
        std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{});
  }

  /**
   * \brief Helper function to forward as a reference wrapper if the type is a reference.
   *
   * \tparam T The type to be forwarded.
   * \param t The value to be forwarded.
   * \return auto The result after forwarding as a reference wrapper if the type is a reference.
   */
  template <typename T>
  auto forwardasReferenceWrapperIfIsReference(T&& t) {
    if constexpr (std::is_lvalue_reference_v<decltype(t)>)
      return std::ref(t);
    else
      return t;
  }

  /**
   * \brief Helper function to make a tuple of values and reference wrappers.
   *
   * \tparam Pars The tuple type.
   * \tparam Tuple The tuple of arguments.
   * \tparam I Index sequence for tuple elements.
   * \param t The tuple of arguments.
   * \param p The tuple of parameters.
   * \return constexpr decltype(auto) The resulting tuple of values and reference wrappers.
   */
  template <class Pars, class Tuple, std::size_t... I>
  constexpr decltype(auto) makeTupleOfValuesAndReferences(Tuple&& t, Pars&& p, std::index_sequence<I...>) {
    return std::make_tuple(
        forwardasReferenceWrapperIfIsReference(applyAndRemoveReferenceWrapper(std::get<I>(t), p.args))...);
  }

  /**
   * \brief Represents a tuple of functions.
   *
   * \tparam Args The argument types.
   */
  template <typename... Args>
  struct Functions
  {
    std::tuple<std::reference_wrapper<std::remove_reference_t<Args>>...> args;
  };

  /**
   * \brief Represents a tuple of parameters.
   *
   * \tparam Args The argument types.
   */
  template <typename... Args>
  struct Parameter
  {
    std::tuple<std::reference_wrapper<std::remove_reference_t<Args>>...> args;
  };

} // namespace Impl

/**
 * \brief Creates a Parameter object.
 *
 * \tparam Args The argument types.
 * \param args The tuple of arguments.
 * \return auto The Parameter object.
 */
template <typename... Args>
auto parameter(Args&&... args) {
  return Impl::Parameter<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

/**
 * \brief Creates a Functions object.
 *
 * \tparam Args The argument types.
 * \param args The tuple of arguments.
 * \return auto The Functions object.
 */
template <typename... Args>
auto functions(Args&&... args) {
  return Impl::Functions<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

/**
 * \brief Initializes the results for functions and parameters.
 *
 * \tparam DerivativeArgs The types of derivative arguments.
 * \tparam ParameterArgs The types of parameter arguments.
 * \param derivativesFunctions The Functions object for derivative arguments.
 * \param parameter The Parameter object for parameter arguments.
 * \return auto The initialized results.
 */
template <typename... DerivativeArgs, typename... ParameterArgs>
auto initResults(const Impl::Functions<DerivativeArgs...>& derivativesFunctions,
                 const Impl::Parameter<ParameterArgs...>& parameter) {
  return Impl::makeTupleOfValuesAndReferences(
      derivativesFunctions.args, parameter,
      std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<decltype(derivativesFunctions.args)>>>{});
}

/**
 * \brief Represents a NonLinearOperator class for handling nonlinear operators.
 * \ingroup utils
 * \tparam TypeListOne The type list for the first set of functions.
 * \tparam TypeListTwo The type list for the second set of functions.
 */
template <typename TypeListOne, typename TypeListTwo>
class NonLinearOperator
{
public:
  NonLinearOperator([[maybe_unused]] const TypeListOne& derivativesFunctions,
                    [[maybe_unused]] const TypeListTwo& args) {
    static_assert(!sizeof(TypeListOne),
                  "This type should not be instantiated. check that your arguments satisfies the template below");
  }
};

/**
 * \brief NonLinearOperator is a class taking linear algebra function and their arguments.
 * The fcuntion are assumed to be derivatvies of each other w.r.t. the first parameter
 *
 * \tparam DerivativeArgs The types of derivative arguments.
 * \tparam ParameterArgs The types of parameter arguments.
 */
template <typename... DerivativeArgs, typename... ParameterArgs>
class NonLinearOperator<Impl::Functions<DerivativeArgs...>, Impl::Parameter<ParameterArgs...>>
{
public:
  using FunctionReturnValues =
      std::tuple<Ikarus::traits::ReturnType<DerivativeArgs, ParameterArgs&...>...>; ///< Function return values
  using ParameterValues = std::tuple<ParameterArgs...>;                             ///< Types of the parameters

  /**
   * \brief Alias for the return type of a function.
   *
   * \tparam n Index of the function.
   */
  template <int n>
  using FunctionReturnType = std::tuple_element_t<n, FunctionReturnValues>;

  /**
   * \brief Alias for the parameter type.
   *
   * \tparam n Index of the parameter.
   */
  template <int n>
  using ParameterValue = std::remove_cvref_t<std::tuple_element_t<n, ParameterValues>>;

  using ValueType =
      std::remove_cvref_t<std::tuple_element_t<0, FunctionReturnValues>>; ///< Return value of the first function
  using DerivativeType =
      std::remove_cvref_t<std::tuple_element_t<1, FunctionReturnValues>>; ///< Return value of the second function

  /**
   * \brief Constructor for NonLinearOperator.
   *
   * \param derivativesFunctions The Functions object for derivative arguments.
   * \param parameterI The Parameter object for parameter arguments.
   */
  explicit NonLinearOperator(const Impl::Functions<DerivativeArgs...>& derivativesFunctions,
                             const Impl::Parameter<ParameterArgs...>& parameterI)
      : derivatives_{derivativesFunctions.args},
        args_{parameterI.args},
        derivativesEvaluated_(initResults(derivativesFunctions, parameterI)) {}

  /**
   * \brief Updates all functions.
   *
   * This function is usually called if the parameters change.
   */
  void updateAll() {
    Dune::Hybrid::forEach(
        Dune::Hybrid::integralRange(Dune::index_constant<sizeof...(DerivativeArgs)>()), [&](const auto i) {
          std::get<i>(derivativesEvaluated_) = Impl::applyAndRemoveReferenceWrapper(std::get<i>(derivatives_), args_);
        });
  }

  /**
   * \brief Updates the n-th function.
   *
   * \tparam n Index of the function to update.
   */
  template <int n>
  void update() {
    std::get<n>(derivativesEvaluated_) = Impl::applyAndRemoveReferenceWrapper(std::get<n>(derivatives_), args_);
  }

  /**
   * \brief Returns the value of the zeroth function.
   *
   * This corresponds to the energy value.
   *
   * \return auto& Reference to the zeroth function value.
   */
  auto& value()
  requires(sizeof...(DerivativeArgs) > 0)
  {
    return nthDerivative<0>();
  }

  /**
   * \brief Returns the derivative value.
   *
   * This corresponds to the gradient of an energy.
   *
   * \return auto& Reference to the derivative function value.
   */
  auto& derivative()
  requires(sizeof...(DerivativeArgs) > 1)
  {
    return nthDerivative<1>();
  }

  /**
   * \brief Returns the second derivative value.
   *
   * This corresponds to the Hessian of an energy.
   *
   * \return auto& Reference to the second derivative function value.
   */
  auto& secondDerivative()
  requires(sizeof...(DerivativeArgs) > 2)
  {
    return nthDerivative<2>();
  }

  /**
   * \brief Returns the n-th derivative value.
   *
   * \tparam n Index of the derivative to return.
   * \return auto& Reference to the n-th derivative function value.
   */
  template <int n>
  auto& nthDerivative()
  requires(sizeof...(DerivativeArgs) > n)
  {
    if constexpr (requires { std::get<n>(derivativesEvaluated_).get(); })
      return std::get<n>(derivativesEvaluated_).get();
    else
      return std::get<n>(derivativesEvaluated_);
  }

  /**
   * \brief Returns the last parameter value.
   *
   * \return auto& Reference to the last parameter value.
   */
  auto& lastParameter() { return nthParameter<sizeof...(ParameterArgs) - 1>(); }
  /**
   * \brief Returns the first parameter value.
   *
   * \return auto& Reference to the first parameter value.
   */
  auto& firstParameter()
  requires(sizeof...(ParameterArgs) > 0)
  {
    return nthParameter<0>();
  }
  /**
   * \brief Returns the second parameter value.
   *
   * \return auto& Reference to the second parameter value.
   */
  auto& secondParameter()
  requires(sizeof...(ParameterArgs) > 1)
  {
    return nthParameter<1>();
  }
  /**
   * \brief Returns the n-th parameter value.
   *
   * \tparam n Index of the parameter to return.
   * \return auto& Reference to the n-th parameter value.
   */
  template <int n>
  auto& nthParameter()
  requires(sizeof...(ParameterArgs) >= n)
  {
    return std::get<n>(args_).get();
  }

  /**
   * \brief Returns a new NonLinearOperator from the given indices.
   *
   * \tparam Derivatives Indices of the functions to include.
   * \return auto The new NonLinearOperator.
   */
  template <int... Derivatives>
  auto subOperator() {
    return Ikarus::NonLinearOperator(functions(std::get<Derivatives>(derivatives_)...),
                                     Impl::applyAndRemoveReferenceWrapper(parameter<ParameterArgs...>, args_));
  }

private:
  using FunctionReturnValuesWrapper = std::tuple<std::conditional_t<
      std::is_reference_v<Ikarus::traits::ReturnType<DerivativeArgs, ParameterArgs&...>>,
      std::reference_wrapper<std::remove_reference_t<Ikarus::traits::ReturnType<DerivativeArgs, ParameterArgs&...>>>,
      std::remove_cvref_t<Ikarus::traits::ReturnType<DerivativeArgs, ParameterArgs&...>>>...>;
  std::tuple<std::conditional_t<std::is_reference_v<DerivativeArgs>,
                                std::reference_wrapper<std::remove_reference_t<DerivativeArgs>>,
                                std::remove_reference_t<DerivativeArgs>>...>
      derivatives_;
  std::tuple<std::conditional_t<std::is_reference_v<ParameterArgs>,
                                std::reference_wrapper<std::remove_reference_t<ParameterArgs>>,
                                std::remove_reference_t<ParameterArgs>>...>
      args_;
  FunctionReturnValuesWrapper derivativesEvaluated_{};
};
} // namespace Ikarus
