//
// Created by Alex on 27.05.2021.
//

#pragma once
#include <concepts>
#include <tuple>
#include <type_traits>
namespace Ikarus::Std {

template <typename> struct is_tuple: std::false_type {};

template <typename ...T> struct is_tuple<std::tuple<T...>>: std::true_type {};


template <class Tuple,class Type> requires is_tuple<Tuple>::value
consteval int countType() {
  int count = 0;
  Dune::Hybrid::forEach(
      Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<Tuple>>()), [&](auto i) {
        using currentType = std::remove_cvref_t<std::tuple_element_t<i, Tuple>>;
        if constexpr (std::is_same_v<currentType, Type>) ++count;
      });
  return count;
}


template <typename T, typename Tuple>
struct hasType;

template <typename T>
struct hasType<T, std::tuple<>> : std::false_type {};

template <typename T, typename U, typename... Ts>
struct hasType<T, std::tuple<U, Ts...>> : hasType<T, std::tuple<Ts...>> {};

template <typename T, typename... Ts>
struct hasType<T, std::tuple<T, Ts...>> : std::true_type {};


template<template<typename...> class, typename...>
struct isInstantiation : public std::false_type {};

template<template<typename...> class U, typename... T>
struct isInstantiation<U, U<T...>> : public std::true_type {};

}  // namespace Ikarus::utils