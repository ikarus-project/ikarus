// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once
// https://www.scs.stanford.edu/~dm/blog/va-opt.html

#define PARENS ()
#define ENUM_CASE(name) \
  case name:            \
    return #name;
// Limits to 256 items
#define EXPAND(arg) EXPAND1(EXPAND1(EXPAND1(EXPAND1(arg))))
#define EXPAND1(arg) EXPAND2(EXPAND2(EXPAND2(EXPAND2(arg))))
#define EXPAND2(arg) EXPAND3(EXPAND3(EXPAND3(EXPAND3(arg))))
#define EXPAND3(arg) EXPAND4(EXPAND4(EXPAND4(EXPAND4(arg))))
#define EXPAND4(arg) arg

#define FOR_EACH(macro, ...) __VA_OPT__(EXPAND(FOR_EACH_HELPER(macro, __VA_ARGS__)))
#define FOR_EACH_HELPER(macro, a1, ...) macro(a1) __VA_OPT__(FOR_EACH_AGAIN PARENS(macro, __VA_ARGS__))
#define FOR_EACH_AGAIN() FOR_EACH_HELPER

#define ENUM_CASE(name) \
  case name:            \
    return #name;

#define MAKE_ENUM(type, ...)                   \
  enum class type { BEGIN, __VA_ARGS__, END }; \
  constexpr std::string toString(type _e) {    \
    using enum type;                           \
    switch (_e) {                              \
      FOR_EACH(ENUM_CASE, __VA_ARGS__)         \
      ENUM_CASE(BEGIN)                         \
      ENUM_CASE(END)                           \
    }                                          \
    __builtin_unreachable();                   \
  }

#include <dune/common/exceptions.hh>
namespace Ikarus {
  template <typename MessageType>
  MessageType& increment(MessageType& e) {
    if (e == MessageType::END) {
      DUNE_THROW(Dune::RangeError, "for MessageType& operator ++ (MessageType&)");
    }
    e = MessageType(static_cast<typename std::underlying_type<MessageType>::type>(e) + 1);
    return e;
  }
}  // namespace Ikarus
