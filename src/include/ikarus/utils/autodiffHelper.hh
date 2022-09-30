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
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace Ikarus {

  template <typename Fun, typename... Vars, typename... Args, typename U, typename G, typename H>
  void hessianN(const Fun &f, const autodiff::Wrt<Vars...> &wrt, const autodiff::At<Args...> &at, U &u,
                std::array<G, U::RowsAtCompileTime> &g, std::array<H, U::RowsAtCompileTime> &h) {
    static_assert(sizeof...(Vars) >= 1);
    static_assert(sizeof...(Args) >= 1);

    auto fEntry = [&](auto &I) { return [&](const auto &) { return std::apply(f, at.args)[I]; }; };
    for (int i = 0; i < U::RowsAtCompileTime; ++i)
      hessian(fEntry(i), wrt, at, u[i], g[i], h[i]);
  }
}  // namespace Ikarus