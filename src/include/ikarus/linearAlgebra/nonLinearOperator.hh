//
// Created by Alex on 21.07.2021.
//

#pragma once
#include <dune/common/hybridutilities.hh>

#include <ikarus/utils/utils/traits.hh>

template <typename Fun, typename... Args>
using ReturnType = std::invoke_result_t<Fun, Args...>;

namespace Impl {
  template <class F, class Tuple, std::size_t... I>
  constexpr decltype(auto) applyAndRemoveRefererenceWrapper(F&& f, Tuple&& t, std::index_sequence<I...>) {
    return std::invoke(
        std::forward<F>(f),
        std::get<I>(std::forward<Tuple>(t)).get()...);  //.get gets the impl type of std::referenceWrapper
  }
}  // namespace Impl

template <class F, class Tuple>
constexpr decltype(auto) applyAndRemoveReferenceWrapper(F&& f, Tuple&& t) {
  return Impl::applyAndRemoveRefererenceWrapper(
      std::forward<F>(f), std::forward<Tuple>(t),
      std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{});
}

template <typename... Args>
struct Parameter {
  std::tuple<std::reference_wrapper<std::remove_cvref_t<Args>>...> args;
};

template <typename... Args>
auto parameter(Args&&... args) {
  return Parameter<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

template <typename... Args>
struct LinearAlgebraFunctions {
  std::tuple<std::reference_wrapper<std::remove_cvref_t<Args>>...> args;
};

template <typename... Args>
auto linearAlgebraFunctions(Args&&... args) {
  return LinearAlgebraFunctions<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

namespace Ikarus {
  template <typename T>
  auto forwardasReferenceWrapperIfIsReference(T&& t) {
    if constexpr (std::is_lvalue_reference_v<decltype(t)>)
      return std::ref(t);
    else
      return t;
  }

  template <class Pars, class Tuple, std::size_t... I>
  constexpr decltype(auto) testF(Tuple&& t, Pars&& p, std::index_sequence<I...>) {
    return std::make_tuple(
        forwardasReferenceWrapperIfIsReference(applyAndRemoveReferenceWrapper(std::get<I>(t), p.args))...);
  }
  template <typename... DerivativeArgs, typename... ParameterArgs>
  auto initResults(const LinearAlgebraFunctions<DerivativeArgs...>& derivativesFunctions,
                   const Parameter<ParameterArgs...>& parameterI) {
    return testF(
        derivativesFunctions.args, parameterI,
        std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<decltype(derivativesFunctions.args)>>>{});
  }

  template <typename TypeListOne, typename TypeListTwo>
  class NonLinearOperator {
  public:
    NonLinearOperator([[maybe_unused]] const TypeListOne& derivativesFunctions,
                      [[maybe_unused]] const TypeListTwo& args) {
      static_assert(!sizeof(TypeListOne),
                    "This type should not be instantiated. check that your arguments satisfies the template below");
    }
  };

  template <typename... DerivativeArgs, typename... ParameterArgs>
  class NonLinearOperator<LinearAlgebraFunctions<DerivativeArgs...>, Parameter<ParameterArgs...>> {
  public:
    explicit NonLinearOperator(const LinearAlgebraFunctions<DerivativeArgs...>& derivativesFunctions,
                               const Parameter<ParameterArgs...>& parameterI)
        : derivatives_{derivativesFunctions.args},
          args_{parameterI.args},
          derivativesEvaluated_(initResults(derivativesFunctions, parameterI)) {
      updateAll();
    }

    void updateAll() {
      Dune::Hybrid::forEach(
          Dune::Hybrid::integralRange(Dune::index_constant<sizeof...(DerivativeArgs)>()), [&](const auto i) {
            std::get<i>(derivativesEvaluated_) = applyAndRemoveReferenceWrapper(std::get<i>(derivatives_), args_);
          });
    }
    template <int n>
    void update() {
      std::get<n>(derivativesEvaluated_) = applyAndRemoveReferenceWrapper(std::get<n>(derivatives_), args_);
    }

  private:
    using FunctionReturnValues = std::tuple<ReturnType<DerivativeArgs, ParameterArgs&...>...>;
    using ParameterValues      = std::tuple<ParameterArgs...>;

  public:
    template <int n>
    using FunctionReturnType = std::tuple_element_t<n, FunctionReturnValues>;

    template <int n>
    using Parameter = std::remove_cvref_t<std::tuple_element_t<n, ParameterValues>>;

    auto& value() requires(sizeof...(DerivativeArgs) > 0) {
      if constexpr (requires { std::get<0>(derivativesEvaluated_).get(); })
        return std::get<0>(derivativesEvaluated_).get();
      else
        return std::get<0>(derivativesEvaluated_);
    }
    auto& derivative() & requires(sizeof...(DerivativeArgs) > 1) {
      if constexpr (requires { std::get<1>(derivativesEvaluated_).get(); })
        return std::get<1>(derivativesEvaluated_).get();
      else
        return std::get<1>(derivativesEvaluated_);
    }
    auto& secondDerivative() requires(sizeof...(DerivativeArgs) > 2) {
      if constexpr (requires { std::get<2>(derivativesEvaluated_).get(); })
        return std::get<2>(derivativesEvaluated_).get();
      else
        return std::get<2>(derivativesEvaluated_);
    }
    template <int n>
    auto& nthDerivative() requires(sizeof...(DerivativeArgs) > n) {
      if constexpr (requires { std::get<n>(derivativesEvaluated_).get(); })
        return std::get<n>(derivativesEvaluated_).get();
      else
        return std::get<n>(derivativesEvaluated_);
    }

    template <int... Derivatives>
    auto subOperator() {
      return Ikarus::NonLinearOperator(linearAlgebraFunctions(std::get<Derivatives>(derivatives_)...),
                                       applyAndRemoveReferenceWrapper(parameter<ParameterArgs...>, args_));
    }

    template <int n>
    auto& nthParameter() requires(sizeof...(ParameterArgs) >= n) {
      return std::get<n>(args_).get();
    }

    auto& lastParameter() { return std::get<sizeof...(ParameterArgs) - 1>(args_).get(); }
    auto& firstParameter() { return std::get<0>(args_).get(); }
    auto& secondParameter() requires(sizeof...(ParameterArgs) > 1) { return std::get<1>(args_).get(); }

    using ValueType      = std::remove_cvref_t<std::tuple_element_t<0, FunctionReturnValues>>;
    using DerivativeType = std::remove_cvref_t<std::tuple_element_t<1, FunctionReturnValues>>;

  private:
    using FunctionReturnValuesWrapper = std::tuple<
        std::conditional_t<std::is_reference_v<ReturnType<DerivativeArgs, ParameterArgs&...>>,
                           std::reference_wrapper<std::remove_cvref_t<ReturnType<DerivativeArgs, ParameterArgs&...>>>,
                           std::remove_cvref_t<ReturnType<DerivativeArgs, ParameterArgs&...>>>...>;
    std::tuple<std::conditional_t<std::is_reference_v<DerivativeArgs>,
                                  std::reference_wrapper<std::remove_cvref_t<DerivativeArgs>>,
                                  std::remove_cvref_t<DerivativeArgs>>...>
        derivatives_;
    std::tuple<std::conditional_t<std::is_reference_v<ParameterArgs>,
                                  std::reference_wrapper<std::remove_cvref_t<ParameterArgs>>,
                                  std::remove_cvref_t<ParameterArgs>>...>
        args_;
    FunctionReturnValuesWrapper derivativesEvaluated_{};
  };
}  // namespace Ikarus