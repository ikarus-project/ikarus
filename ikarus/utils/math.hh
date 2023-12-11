// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once
#include <concepts>

namespace Ikarus {
  namespace Impl {
    template <typename T>
    constexpr T sqrt_helper(T x, T lo, T hi) {
      if (lo == hi) return lo;

      const T mid = (lo + hi + 1) / 2;

      if (x / mid < mid)
        return sqrt_helper<T>(x, lo, mid - 1);
      else
        return sqrt_helper(x, mid, hi);
    }

  }  // namespace Impl

  // compile time sqrt for integer types
  template <typename T>
  requires std::integral<T>
  constexpr T ct_sqrt(T x) { return Impl::sqrt_helper<T>(x, 0, x / 2 + 1); }
}  // namespace Ikarus
