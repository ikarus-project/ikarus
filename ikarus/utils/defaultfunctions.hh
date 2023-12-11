// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

namespace Ikarus {
  struct SolverDefault {
    template <typename A, typename B>
    constexpr auto operator()(A&& a, B&& b) const {
      return a / b;
    }
  };

  struct UpdateDefault {
    template <typename A, typename B>
    constexpr void operator()(A&& a, B&& b) const {
      a += b;
    }
  };

  struct LoadDefault {};

}  // namespace Ikarus
