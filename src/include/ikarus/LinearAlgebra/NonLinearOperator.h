//
// Created by Alex on 21.07.2021.
//

#pragma once
#include <ikarus/utils/utils/traits.h>
#include <dune/common/hybridutilities.hh>

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
struct Derivatives {
  std::tuple<std::reference_wrapper<std::remove_cvref_t<Args>>...>  args;
};

template <typename... Args>
auto derivatives(Args&&... args) {
  return Derivatives<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

namespace Ikarus {
  template <typename TypeListOne, typename TypeListTwo, typename TypeListThree>
  class NonLinearOperator {
  public:
    NonLinearOperator(const TypeListOne& valueFunc,
                      const TypeListTwo& derivativesFunctions,
                      const TypeListThree& args){}

  };

  template <typename ValueFuncType, typename... DerivativeArgs, typename... ParameterArgs>
  class NonLinearOperator<ValueFuncType, Derivatives<DerivativeArgs...>, Parameter<ParameterArgs...>> {
  public:

     NonLinearOperator(const ValueFuncType& valueFunc,
                               const Derivatives<DerivativeArgs...>& derivativesFunctions,
                               const Parameter<ParameterArgs...>& parameter)
    : valueFunction_{valueFunc},
      derivatives_{derivativesFunctions.args} ,
       args_{parameter.args}{
           updateAll();
    }

    void updateAll() {
      value_                             = std::apply(valueFunction_, args_);
      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<sizeof...(DerivativeArgs)>()), [&](const auto i) {
        std::get<i>(derivativesEvaluated_) = std::apply(std::get<i>(derivatives_), args_);
      });
    }
    template<int n>
    void update()
    {
      if constexpr (n==0)
        value_                             = std::apply(valueFunction_, args_);
      else
        std::get<n-1>(derivativesEvaluated_) = std::apply(std::get<n-1>(derivatives_), args_);
    }
    using ValueReturnType = ReturnType<ValueFuncType, ParameterArgs&...>;

    ValueReturnType& value() { return value_; }
    auto& derivative() requires(sizeof...(DerivativeArgs) > 0) { return std::get<0>(derivativesEvaluated_); }
    auto& secondDerivative() requires(sizeof...(DerivativeArgs) > 1) { return std::get<1>(derivativesEvaluated_); }
    template <int n>
    auto& nthDerivative() requires(sizeof...(DerivativeArgs) > n - 1) {
      return std::get<n - 1>(derivativesEvaluated_);
    }

  private :
      ValueFuncType valueFunction_;
    std::tuple<std::reference_wrapper<std::remove_cvref_t<DerivativeArgs>>...>  derivatives_;
    std::tuple<std::reference_wrapper<std::remove_cvref_t<ParameterArgs>>...>  args_;
    ValueReturnType value_;
    std::tuple<ReturnType<DerivativeArgs, ParameterArgs&...>...> derivativesEvaluated_;
  };
}  // namespace Ikarus