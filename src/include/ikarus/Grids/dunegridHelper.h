//
// Created by Alex on 26.01.2022.
//

#pragma once

#include <ranges>

#include <dune/geometry/dimension.hh>
#include <dune/grid/common/partitionset.hh>
namespace Dune {

  template <typename GV>
    requires requires(GV gv) { entities(gv, Dune::Dim<3>(), Dune::Partitions::all); }
  auto volumes(GV& gv) {
    return entities(gv, Dune::Dim<3>(), Dune::Partitions::all);
  }

  template <typename GV>
    requires requires(GV gv) { entities(gv, Dune::Dim<3>(), Dune::Partitions::all); }
  auto volumes(GV const& gv) {
    return entities(gv, Dune::Dim<3>(), Dune::Partitions::all);
  }

  template <typename GV>
    requires requires(GV gv) { entities(gv, Dune::Dim<2>(), Dune::Partitions::all); }
  auto surfaces(GV& gv) {
    return entities(gv, Dune::Dim<2>(), Dune::Partitions::all);
  }

  template <typename GV>
    requires requires(GV gv) { entities(gv, Dune::Dim<2>(), Dune::Partitions::all); }
  auto surfaces(GV const& gv) {
    return entities(gv, Dune::Dim<2>(), Dune::Partitions::all);
  }

  template <typename GV>
    requires requires(GV gv) { entities(gv, Dune::Dim<1>(), Dune::Partitions::all); }
  auto edges(GV& gv) {
    return entities(gv, Dune::Dim<1>(), Dune::Partitions::all);
  }

  template <typename GV>
    requires requires(GV gv) { entities(gv, Dune::Dim<1>(), Dune::Partitions::all); }
  auto edges(GV const& gv) {
    return entities(gv, Dune::Dim<1>(), Dune::Partitions::all);
  }

  template <typename GE>
    requires requires(GE ge, int i) { ge.template subEntity<GE::dimension - 1>(i); }
  auto edges(GE const& ge) {
    return std::views::iota(0U, ge.subEntities(GE::dimension - 1))
           | std::views::transform([&](auto&& i) { return ge.template subEntity<GE::dimension - 1>(i); });
  }

  template <typename GE>
    requires requires(GE ge, int i) { ge.subEntity<GE::dimension - 1>(i); }
  auto edges(GE& ge) {
    return std::views::iota(0U, ge.subEntities(GE::dimension - 1))
           | std::views::transform([&](auto&& i) { return ge.template subEntity<GE::dimension - 1>(i); });
  }

  inline Dune::GeometryType duneType(Dune::GeometryType duneType) { return duneType; }

}  // namespace Dune
