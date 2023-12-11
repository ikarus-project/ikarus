// Original File: https://gitlab.mn.tu-dresden.de/amdis/amdis-core/-/blob/master/amdis/functions/FlatPreBasis.hpp
// SPDX-FileCopyrightText: 2023 Copyright (c)  AMDiS
// SPDX-License-Identifier: MIT
// Modifications:
// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <cstddef>
#include <utility>

#include <dune/common/indices.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

namespace Ikarus {
  /// \brief Transform a PreBasis into one with flat index-merging strategy
  /**
   * This utility takes a pre-basis and converts recursively all index-merging strategies
   * into their flat analog, i.e. `BlockedInterleaved` is converted into `FlatInterleaved`
   * and `BlockedLexicographic` is transformed into `FlatLexicographic`.
   *
   * This type-trait needs to be specialized for all PreBasis types that need special handling,
   * like PowerPreBasis or CompositePreBasis.
   *
   * \relates flatPreBasis()
   * \relates FlatPreBasis_t
   **/
  template <class PreBasis>
  struct FlatPreBasis {
    using type = PreBasis;

    /// Try to construct the pre-basis using a `gridView`.
    template <class PB>
    static type create(PB const& preBasis) {
      return {preBasis.gridView()};
    }

    /// Do not transform the preBasis if already flat
    static PreBasis const& create(PreBasis const& preBasis) { return preBasis; }
  };

  /// \brief Type alias for flatted PreBasis
  template <class PreBasis>
  using FlatPreBasis_t = typename FlatPreBasis<PreBasis>::type;

  /// \brief Define the flat index-merging strategy for a given strategy `IMS`
  template <class IMS>
  struct FlatIndexMergingStrategy {
    using type = IMS;
  };

  // specialization for BlockedInterleaved
  template <>
  struct FlatIndexMergingStrategy<Dune::Functions::BasisFactory::BlockedInterleaved> {
    using type = Dune::Functions::BasisFactory::FlatInterleaved;
  };

  // specialization for BlockedLexicographic
  template <>
  struct FlatIndexMergingStrategy<Dune::Functions::BasisFactory::BlockedLexicographic> {
    using type = Dune::Functions::BasisFactory::FlatLexicographic;
  };

  // specialization for composite bases
  template <class IMS, class... SPB>
  struct FlatPreBasis<Dune::Functions::CompositePreBasis<IMS, SPB...>> {
    using FIMS = typename FlatIndexMergingStrategy<IMS>::type;
    using type = Dune::Functions::CompositePreBasis<FIMS, FlatPreBasis_t<SPB>...>;

    template <class PreBasis>
    static type create(PreBasis const& preBasis) {
      return create(preBasis, std::index_sequence_for<SPB...>{});
    }

    template <class PreBasis, std::size_t... I>
    static type create(PreBasis const& preBasis, std::index_sequence<I...>) {
      return {FlatPreBasis<SPB>::create(preBasis.subPreBasis(Dune::index_constant<I>{}))...};
    }
  };

  // specialization for power bases
  template <class IMS, class SPB, std::size_t C>
  struct FlatPreBasis<Dune::Functions::PowerPreBasis<IMS, SPB, C>> {
    using FIMS = typename FlatIndexMergingStrategy<IMS>::type;
    using type = Dune::Functions::PowerPreBasis<FIMS, FlatPreBasis_t<SPB>, C>;

    template <class PreBasis>
    static type create(PreBasis const& preBasis) {
      return {FlatPreBasis<SPB>::create(preBasis.subPreBasis())};
    }
  };

  /// \brief Generator function for a flatted PreBasis
  template <class PreBasis>
  decltype(auto) flatPreBasis(PreBasis const& preBasis) {
    return FlatPreBasis<PreBasis>::create(preBasis);
  }

}  // end namespace Ikarus
