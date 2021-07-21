//
// Created by Alex on 21.07.2021.
//

#pragma once

template<typename Fun, typename... Args>
using ReturnType = std::invoke_result_t<Fun, Args...>;


template <typename DerivFuncType, typename ValueFuncType, typename... Args>
class NonLinearOperator {
public:
  NonLinearOperator(const  ValueFuncType& valueFunc,const  DerivFuncType& jacobianFunc, Args&... args)
      : valueFunction_{valueFunc}, jacobianFunction_{jacobianFunc}, args_{std::forward_as_tuple(std::forward<Args&>(args)...)}
  {
    updateAll();
  }

  void updateAll() {
    value_ = std::apply(valueFunction_,args_);
    deriv_ = std::apply(jacobianFunction_,args_);
  }
  using ValueReturnType = ReturnType<ValueFuncType,Args&...>;
  using DerivativeReturnType = ReturnType<DerivFuncType,Args&...>;

  ValueReturnType& value() { return value_; }
  DerivativeReturnType& derivative() { return deriv_; }

private:

  ValueFuncType valueFunction_;
  DerivFuncType jacobianFunction_;
  ValueReturnType value_;
  DerivativeReturnType deriv_;  // gradient or jacobian
  std::tuple<std::reference_wrapper<Args>...> args_;
};
