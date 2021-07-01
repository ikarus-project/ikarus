//
// Created by Alex on 25.05.2021.
//

#pragma once

// This file contains stl-like algorithms
#include <ranges>
namespace Ikarus::stl{
void makeUniqueAndSort(std::ranges::random_access_range auto& varVec) {
  sort(varVec.begin(), varVec.end());
  varVec.erase(std::unique(varVec.begin(), varVec.end()), varVec.end());
}

template <typename Value>
auto appendUnique(std::ranges::random_access_range auto& c, Value&& v) {
  static_assert(std::is_same_v<typename decltype(begin(c))::value_type, std::remove_reference_t<decltype(v)>>);
  const auto it = find(begin(c), end(c), v);
  size_t index = std::distance(begin(c), it);
  if (it == end(c)) c.push_back(std::forward<Value>(v));

  assert(!c.empty());
  return index;
}

template <class Container>  // TODO: create concept for this
void printContent(Container& varVec) {
  std::ranges::for_each(varVec, [](auto& var) { std::cout << var << '\n'; });
}

template <class Container>
auto transformValueRangeToPointerRange(Container& varVec)
{
  auto transformValueToPointer = [](auto&& obj){ return &obj;};
  return  (std::ranges::subrange(varVec.begin(), varVec.end()) | std::views::transform(transformValueToPointer));
}

}
