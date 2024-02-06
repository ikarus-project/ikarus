// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file traits.hh
 * \brief Contains stl-like type traits
 */

#pragma once
#include <optional>
#include <tuple>
#include <type_traits>

#include <dune/common/hybridutilities.hh>
namespace Ikarus::traits {

/**
 * \brief Concept to check if a type is a pointer or nullptr_t.
 * \ingroup traits
 * \tparam T Type to check.
 */
template <typename T>
concept Pointer = std::is_pointer_v<T> || std::is_same_v<T, std::nullptr_t>;

#ifndef DOXYGEN
template <typename>
struct is_tuple : std::false_type
{
};
#endif
/**
 * \brief Type trait to check if a type is an instantiation of std::tuple.
 * \ingroup traits
 *
 */
template <typename... T>
struct is_tuple<std::tuple<T...>> : std::true_type
{
};

/**
 * \brief Metafunction to count the occurrences of a specific type in a tuple.
 * \ingroup traits
 * \tparam Tuple Type of the tuple.
 * \tparam Type Type to count in the tuple.
 * \return int Number of occurrences of the specified type in the tuple.
 */
template <class Tuple, class Type>
requires is_tuple<Tuple>::value
consteval int countType() {
  int count = 0;
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<Tuple>>()), [&](auto i) {
    using currentType = std::remove_cvref_t<std::tuple_element_t<i, Tuple>>;
    if constexpr (std::is_same_v<currentType, Type>)
      ++count;
  });
  return count;
}

/**
 * \brief Type trait to obtain the return type of a callable type when given specific arguments.
 *
 * \ingroup traits
 *
 * \tparam Fun Callable type.
 * \tparam Args Argument types.
 */
template <typename Fun, typename... Args>
using ReturnType = std::invoke_result_t<Fun, Args...>;

/**
 * \brief Type trait to check if a specified type is present in a tuple.
 *
 * \ingroup traits
 *
 * \tparam T Type to check.
 * \tparam Tuple Tuple to search within.
 */
template <typename T, typename Tuple>
struct hasType : std::false_type
{
};

#ifndef DOXYGEN
template <typename T>
struct hasType<T, std::tuple<>> : std::false_type
{
};

template <typename T>
struct hasType<T, T> : std::true_type
{
};

/**
 * \brief Recursive template specialization of hasType trait for tuples.
 * \ingroup traits
 * \tparam T Type to check.
 * \tparam U Current tuple type.
 * \tparam Ts Remaining types in the tuple.
 */
template <typename T, typename U, typename... Ts>
struct hasType<T, std::tuple<U, Ts...>> : hasType<T, std::tuple<Ts...>>
{
};

template <typename T, typename... Ts>
struct hasType<T, std::tuple<T, Ts...>> : std::true_type
{
};
#endif

#ifndef DOXYGEN
template <template <typename...> class, typename...>
struct isSpecialization : std::false_type
{
};
#endif

/**
 * \brief Type trait to check if a class is a specialization of a template.
 *
 * \ingroup traits
 *
 */
template <template <typename...> class U, typename... T>
struct isSpecialization<U, U<T...>> : std::true_type
{
};

#ifndef DOXYGEN
template <template <typename, auto...> class Type, typename>
struct isSpecializationTypeAndNonTypes : std::false_type
{
};

template <template <typename, auto...> class Type, typename T, auto... N>
struct isSpecializationTypeAndNonTypes<Type, Type<T, N...>> : std::true_type
{
};

template <template <auto, typename...> class Type, typename>
struct isSpecializationNonTypeAndTypes : std::false_type
{
};
#endif

/**
 * \brief Type trait to check if a class is a specialization of a template with a non-type parameter and types.
 *
 * \ingroup traits
 *
 * \tparam Type Template class with a non-type parameter and types.
 * \tparam T Non-type parameter.
 * \tparam N Types used to instantiate the template.
 */
template <template <auto, typename...> class Type, auto T, typename... N>
struct isSpecializationNonTypeAndTypes<Type, Type<T, N...>> : std::true_type
{
};

#ifndef DOXYGEN
template <template <typename, auto, typename> class Type, typename>
struct isSpecializationTypeNonTypeAndType : std::false_type
{
};
#endif
/**
 * \brief Type trait to check if a class is a specialization of a template with types and two non-type parameters.
 *
 * \ingroup traits
 *
 * \tparam Type Template class with types and two non-type parameters.
 * \tparam T First type parameter.
 * \tparam M First non-type parameter.
 * \tparam N Second type parameter.
 */
template <template <typename, auto, typename> class Type, typename T, auto M, typename N>
struct isSpecializationTypeNonTypeAndType<Type, Type<T, M, N>> : std::true_type
{
};

#ifndef DOXYGEN
template <template <auto...> class Type, typename>
struct isSpecializationNonTypes : std::false_type
{
};
#endif

/**
 * \brief Type trait to check if a class is a specialization of a template with non-type parameters.
 *
 * \ingroup traits
 *
 * \tparam Type Template class with non-type parameters.
 * \tparam N Non-type parameters.
 */
template <template <auto...> class Type, auto... N>
struct isSpecializationNonTypes<Type, Type<N...>> : std::true_type
{
};

/**
 * \brief Type trait to get the index of a type in a tuple.
 * \ingroup traits
 * \details
 * Usage:
 * ```cpp
 * using foo_t = std::tuple<int, double, float>;
 * static_assert(Index<int, foo_t>::value == 0);
 * static_assert(Index<double, foo_t>::value == 1);
 * static_assert(Index<float, foo_t>::value == 2);
 * static_assert(Index<long double, foo_t>::value == 3);
 * ```
 * If the type is not found, the returned index is the size of the tuple.
 *
 * \tparam T Type to find in the tuple.
 * \tparam Tuple Tuple type.
 */
template <class T, class Tuple>
struct Index;
#ifndef DOXYGEN
template <class T>
struct Index<T, std::tuple<>>
{
  static const std::size_t value = 0;
};

template <class T, class... Types>
struct Index<T, std::tuple<T, Types...>>
{
  static constexpr std::size_t value = 0;
};

template <class T, class U, class... Types>
struct Index<T, std::tuple<U, Types...>>
{
  static const std::size_t value = 1 + Index<T, std::tuple<Types...>>::value;
};
#endif

/**
 * \brief Type trait to rebind the underlying type of containers.
 *
 * \details
 * Specialization for types like std::vector<...> and nested std::vector<std::vector>.
 * ```cpp
 * Rebind<std::vector<int>, double>::other; // --> std::vector<double>
 * Rebind<std::vector<std::vector<int>>, double>::other; // --> std::vector<std::vector<double>>
 * ```
 * Specialization for types like std::array<...,N>.
 * ```cpp
 * Rebind<std::array<int, 5>, double>::other; // --> std::array<double, 5>
 * ```
 *
 * \tparam Container Original container type.
 * \tparam NewType New type to rebind to.
 */
template <class Container, class NewType>
struct Rebind;

#ifndef DOXYGEN
/*
 * Specialization for types like std::vector<...> and nested std::vector<std::vector>
 */
template <class OldType, class... Args, template <class...> class Container, class NewType>
struct Rebind<Container<OldType, Args...>, NewType>
{
  using other = Container<NewType, typename Rebind<Args, NewType>::other...>;
};

/*
 * Specialization for types like std::array<...,N>
 */
template <class OldType, std::size_t N, template <class, std::size_t> class Container, class NewType>
struct Rebind<Container<OldType, N>, NewType>
{
  using other = Container<NewType, N>;
};

#endif

/**
 * \ingroup traits
 * \brief Type trait for extracting information about functions.
 *
 * \details
 * This trait provides information about the return type, argument types, and the number of arguments of a function.
 *
 * \tparam T Type of the function.
 */
template <typename T, typename = void>
struct FunctionTraits;

#ifndef DOXYGEN
/**
 * \brief Specialization for general functions
 */
template <typename R, typename... Args>
struct FunctionTraits<R (*)(Args...)>
{
  using return_type = R;
  template <int i>
  using args_type                        = typename std::tuple_element<i, std::tuple<Args...>>::type;
  static constexpr int numberOfArguments = sizeof...(Args);
};

/**
 * \brief Specialization for const member functions
 */
template <typename R, typename C, typename... Args>
struct FunctionTraits<R (C::*)(Args...) const>
{
  using return_type = R;
  template <int i>
  using args_type                        = typename std::tuple_element<i, std::tuple<Args...>>::type;
  static constexpr int numberOfArguments = sizeof...(Args);
};

/**
 * \brief Specialization for non-const member functions
 */
template <typename R, typename C, typename... Args>
struct FunctionTraits<R (C::*)(Args...)>
{
  using return_type = R;
  template <int i>
  using args_type                        = typename std::tuple_element<i, std::tuple<Args...>>::type;
  static constexpr int numberOfArguments = sizeof...(Args);
};

/**
 * \brief Specialization for lambdas using std::void_t to allow a specialization of the original template
 *        The lambda is forwarded using lambdas operator() to the general function traits
 */
template <typename T>
struct FunctionTraits<T, Dune::void_t<decltype(&T::operator())>> : public FunctionTraits<decltype(&T::operator())>
{
};
#endif

} // namespace Ikarus::traits
