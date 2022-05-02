//
// Created by Alex on 02.05.2022.
//

#pragma once
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace Ikarus {

  template <typename Fun, typename... Vars, typename... Args, typename U, typename G, typename H>
  void hessianN(const Fun &f, const autodiff::Wrt<Vars...> &wrt, const autodiff::At<Args...> &at, U &u,
                std::array<G, U::RowsAtCompileTime> &g, std::array<H, U::RowsAtCompileTime> &h) {
    static_assert(sizeof...(Vars) >= 1);
    static_assert(sizeof...(Args) >= 1);

    auto fEntry = [&](auto& I) {
      return  [&](const auto &x) {
        return std::apply(f, at.args)[I];
      };
    };
    for (int i = 0; i <  U::RowsAtCompileTime; ++i)
      hessian(fEntry(i), wrt, at, u[i], g[i], h[i]);
  }
}  // namespace Ikarus