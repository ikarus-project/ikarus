//
// Created by Alex on 27.05.2021.
//

#pragma once
#include <tuple>
#include <concepts>
#include <type_traits>
namespace Ikarus::utils {

  template <typename T>
  struct is_std_array : std::false_type {};

  template <typename T, std::size_t N>
  struct is_std_array<std::array<T, N>> : std::true_type {};

  template <typename T>
  struct function_traits;

  template <typename R, typename... A>
  struct function_traits<R (*)(A...)> {
    using return_type           = R;
    using class_type            = void;
    using args_type             = std::tuple<A...>;
    static constexpr auto arity = sizeof...(A);
  };

  template <typename R, typename... A>
  struct function_traits<R (*&)(A...)> {
    using return_type           = R;
    using class_type            = void;
    using args_type             = std::tuple<A...>;
    static constexpr auto arity = sizeof...(A);
  };

  template <typename R, typename... A>
  struct function_traits<R(A...)> {
    using return_type           = R;
    using class_type            = void;
    using args_type             = std::tuple<A...>;
    static constexpr auto arity = sizeof...(A);
  };

  template <class F, std::size_t... Is, class T>
  std::function<typename T::result_type(std::tuple_element_t<Is, typename T::arg_tuple>...)> lambda_to_func_impl(
      F&& f, std::index_sequence<Is...>, T) {
    return std::function<typename T::result_type(std::tuple_element_t<Is, typename T::arg_tuple>...)>(
        std::forward<F>(f));
  }
  //
  template <class F>
  auto lambda_to_func(F&& f) {
    using traits = function_traits<F>;
    return lambda_to_func_impl(std::forward<F>(f), std::make_index_sequence<traits::arity>{}, traits{});
  }
}  // namespace Ikarus::utils