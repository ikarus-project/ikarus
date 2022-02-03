//
// Created by lex on 02/02/2022.
//

#pragma once
#include <dune/functions/functionspacebases/basistags.hh>

namespace Ikarus::Concepts {

  template <typename Basis>
  concept FlatInterLeavedBasis = requires {
    std::is_same_v<typename Basis::PreBasis::IndexMergingStrategy, Dune::Functions::BasisFactory::FlatInterleaved>;
  };

  template <typename Basis>
  concept FlatLexicographicBasis = requires {
    std::is_same_v<typename Basis::PreBasis::IndexMergingStrategy, Dune::Functions::BasisFactory::FlatLexicographic>;
  };

  template <typename Basis>
  concept FlatIndexBasis = FlatLexicographicBasis<Basis> or FlatInterLeavedBasis<Basis>;

  template <typename Basis>
  concept BlockedInterLeavedBasis = requires {
    std::is_same_v<typename Basis::PreBasis::IndexMergingStrategy, Dune::Functions::BasisFactory::BlockedInterleaved>;
  };

  template <typename Basis>
  concept BlockedLexicographicBasis = requires {
    std::is_same_v<typename Basis::PreBasis::IndexMergingStrategy, Dune::Functions::BasisFactory::BlockedLexicographic>;
  };

  template <typename Basis>
  concept BlockedIndexBasis = BlockedLexicographicBasis<Basis> or BlockedInterLeavedBasis<Basis>;

  template <typename Basis>
  concept PowerBasis = requires {
    Basis::PreBasis::Node::isPower==true;
  };
}  // namespace Ikarus::Concepts
