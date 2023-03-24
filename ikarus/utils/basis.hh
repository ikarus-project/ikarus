//
// Created by lex on 3/23/23.
//

#ifndef IKARUS_BASIS_HH
#define IKARUS_BASIS_HH

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

#include <ikarus/utils/flatPreBasis.hh>

namespace Ikarus {
  template <typename BlockedBasis_, typename FlatBasis_>
  class Basis {
  public:
    using BlockedBasis = std::remove_cvref_t<BlockedBasis_>;
    using GridView     = typename BlockedBasis::GridView;
    using FlatBasis    = std::remove_cvref_t<FlatBasis_>;

    Basis(BlockedBasis_& p_bb, FlatBasis_& p_fb) : bb{p_bb}, fb{p_fb} {}

    auto& flat() { return fb; }

    auto& basis() { return bb; }

    const auto& flat() const { return fb; }

    const auto& basis() const { return bb; }

    const auto& gridView() const { return bb.gridView(); }

    auto& gridView() { return bb.gridView(); }

  private:
    BlockedBasis bb;
    FlatBasis fb;
  };

  template <typename GridView, typename PreBasis>
  auto makeBasis(GridView& gv, const PreBasis& pb) {
    auto basis     = Dune::Functions::BasisFactory::makeBasis(gv, pb);
    auto flatBasis = Dune::Functions::DefaultGlobalBasis(Ikarus::flatPreBasis(pb(gv)));
    return Basis(basis, flatBasis);
  }
}  // namespace Ikarus

#endif  // IKARUS_BASIS_HH
