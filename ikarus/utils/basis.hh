// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <utility>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

#include <ikarus/utils/flatprebasis.hh>

namespace Ikarus {
  template <typename PreBasis_>
  class Basis {
  public:
    using PreBasis       = PreBasis_;
    using GridView       = typename PreBasis::GridView;
    using UntouchedBasis = decltype(Dune::Functions::DefaultGlobalBasis(std::declval<PreBasis>()));
    using FlatBasis = decltype(Dune::Functions::DefaultGlobalBasis(Ikarus::flatPreBasis(std::declval<PreBasis>())));

    explicit Basis(const PreBasis& pb)
        : bb{Dune::Functions::DefaultGlobalBasis(pb)},
          fb{Dune::Functions::DefaultGlobalBasis(Ikarus::flatPreBasis(pb))} {}

    auto& flat() { return fb; }

    auto& untouched() { return bb; }

    const auto& flat() const { return fb; }

    const auto& untouched() const { return bb; }

    const auto& gridView() const { return bb.gridView(); }

    auto& gridView() { return bb.gridView(); }

  private:
    UntouchedBasis bb;
    FlatBasis fb;
  };

  template <typename GridView, typename PreBasisFactory>
  auto makeBasis(const GridView& gv, const PreBasisFactory& pb) {
    auto preBasis = pb(gv);
    return Basis(preBasis);
  }

  template <typename PreBasis>
  auto makeBasis(const Dune::Functions::DefaultGlobalBasis<PreBasis>& gb) {
    return Basis(gb.preBasis());
  }
}  // namespace Ikarus
