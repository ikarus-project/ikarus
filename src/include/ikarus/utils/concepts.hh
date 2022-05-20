//
// Created by lex on 02/02/2022.
//

#pragma once
#include <dune/functions/functionspacebases/basistags.hh>

#include <Eigen/Core>
namespace Ikarus {

  template <typename Derived>
  auto transpose(const Eigen::EigenBase<Derived>& A);
  namespace Concepts {

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
      std::is_same_v<typename Basis::PreBasis::IndexMergingStrategy,
                     Dune::Functions::BasisFactory::BlockedLexicographic>;
    };

    template <typename DuneLocalBasisImpl>
    concept DuneLocalBasis = requires(DuneLocalBasisImpl& duneLocalBasis) {
      typename DuneLocalBasisImpl::Traits::RangeType;
      typename DuneLocalBasisImpl::Traits::JacobianType;
      DuneLocalBasisImpl::Traits::dimDomain;
      typename DuneLocalBasisImpl::Traits::DomainType;

      typename DuneLocalBasisImpl::Traits::DomainFieldType;
      typename DuneLocalBasisImpl::Traits::RangeFieldType;

      duneLocalBasis.evaluateFunction(std::declval<typename DuneLocalBasisImpl::Traits::DomainType>(),
                                      std::declval<std::vector<typename DuneLocalBasisImpl::Traits::RangeType>&>());
      duneLocalBasis.evaluateJacobian(std::declval<typename DuneLocalBasisImpl::Traits::DomainType>(),
                                      std::declval<std::vector<typename DuneLocalBasisImpl::Traits::JacobianType>&>());
    };

    template <typename Basis>
    concept BlockedIndexBasis = BlockedLexicographicBasis<Basis> or BlockedInterLeavedBasis<Basis>;

    template <typename Basis>
    concept PowerBasis = requires {
      Basis::PreBasis::Node::isPower == true;
    };

    template <typename L, typename R>
    concept MultiplyAble = requires(L x, R y) {
      x* y;
    };

    template <typename L, typename R>
    concept AddAble = requires(L x, R y) {
      x + y;
    };

    template <typename L, typename R>
    concept SubstractAble = requires(L x, R y) {
      x - y;
    };

    template <typename L, typename R>
    concept MultiplyAssignAble = requires(L x, R y) {
      x *= y;
    };

    template <typename L, typename R>
    concept DivideAssignAble = requires(L x, R y) {
      x /= y;
    };

    template <typename L, typename R>
    concept AddAssignAble = requires(L x, R y) {
      x += y;
    };

    template <typename L, typename R>
    concept SubstractAssignAble = requires(L x, R y) {
      x -= y;
    };

    template <typename L, typename R>
    concept DivideAble = requires(L x, R y) {
      x / y;
    };

    template <typename L>
    concept NegateAble = requires(L x) {
      -x;
    };

    template <typename L>
    concept TransposeAble = requires(L x) {
      transpose(x);
    };
  }  // namespace Concepts
}  // namespace Ikarus
