/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#pragma once

// This file contains stl-like algorithms
#include <iosfwd>
#include <ranges>
namespace Ikarus::utils {
  void makeUniqueAndSort(std::ranges::random_access_range auto& varVec) {
    sort(varVec.begin(), varVec.end());
    varVec.erase(std::unique(varVec.begin(), varVec.end()), varVec.end());
  }

  /** If the value is not in the range, it is appended*/
  template <typename Value>
  auto appendUnique(std::ranges::random_access_range auto& c, Value&& v) {
    static_assert(std::is_same_v<typename decltype(begin(c))::value_type, std::remove_reference_t<decltype(v)>>);
    const auto it = find(begin(c), end(c), v);
    size_t index  = std::distance(begin(c), it);
    if (it == end(c)) c.push_back(std::forward<Value>(v));

    return index;
  }

  template <class Container>  // TODO: add concept for this
  void printContent(Container&& varVec, std::ostream& os = std::cout) {
    std::ranges::for_each(varVec, [&os](auto&& var) { os << var << '\n'; });
  }

  template <class Container>
  auto transformValueRangeToPointerRange(Container& cont) {
    auto transformValueToPointer = [](auto&& obj) { return &obj; };
    return (std::ranges::subrange(cont.begin(), cont.end()) | std::views::transform(transformValueToPointer));
  }

  template <class Container>
  auto transformPointerRangeToReferenceRange(Container& cont) {
    auto transformValueToPointer = [](auto&& obj) -> auto& { return *obj; };
    return (std::ranges::subrange(cont.begin(), cont.end()) | std::views::transform(transformValueToPointer));
  }

}  // namespace Ikarus::utils
