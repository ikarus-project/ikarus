//
// Created by Alex on 21.07.2021.
//

#pragma once
#include <dune/common/hybridutilities.hh>

#include <ikarus/utils/utils/traits.h>

template <typename Fun, typename... Args>
using ReturnType = std::invoke_result_t<Fun, Args...>;

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
  template <typename TypeListTwo, typename TypeListThree>
  class NonLinearOperator {
  public:
    NonLinearOperator(const TypeListTwo& derivativesFunctions, const TypeListThree& args) {}
  };

  template <typename... DerivativeArgs, typename... ParameterArgs>
  class NonLinearOperator<LinearAlgebraFunctions<DerivativeArgs...>, Parameter<ParameterArgs...>> {
  public:
    explicit NonLinearOperator(const LinearAlgebraFunctions<DerivativeArgs...>& derivativesFunctions,
                               const Parameter<ParameterArgs...>& parameterI)
        : derivatives_{derivativesFunctions.args}, args_{parameterI.args} {
      updateAll();
    }

    void updateAll() {
      Dune::Hybrid::forEach(
          Dune::Hybrid::integralRange(Dune::index_constant<sizeof...(DerivativeArgs)>()),
          [&](const auto i) { std::get<i>(derivativesEvaluated_) = std::apply(std::get<i>(derivatives_), args_); });
    }
    template <int n>
    void update() {
      std::get<n>(derivativesEvaluated_) = std::apply(std::get<n>(derivatives_), args_);
    }
    auto& value() requires(sizeof...(DerivativeArgs) > 0) { return std::get<0>(derivativesEvaluated_); }
    auto& derivative() requires(sizeof...(DerivativeArgs) > 1) { return std::get<1>(derivativesEvaluated_); }
    auto& secondDerivative() requires(sizeof...(DerivativeArgs) > 2) { return std::get<2>(derivativesEvaluated_); }
    template <int n>
    auto& nthDerivative() requires(sizeof...(DerivativeArgs) > n) {
      return std::get<n>(derivativesEvaluated_);
    }

  private : std::tuple<std::reference_wrapper<std::remove_cvref_t<DerivativeArgs>>...> derivatives_;
    std::tuple<std::reference_wrapper<std::remove_cvref_t<ParameterArgs>>...> args_;
    std::tuple<ReturnType<DerivativeArgs, ParameterArgs&...>...> derivativesEvaluated_;
  };
}  // namespace Ikarus