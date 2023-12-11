// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "pathfollowingfunctions.hh"

#include <concepts>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/functions/functionspacebases/basistags.hh>

#include <ikarus/utils/traits.hh>

namespace Eigen {
  template <typename Derived>
  struct EigenBase;
}

namespace Ikarus {
  template <auto stressIndexPair, typename MaterialImpl>
  struct VanishingStress;

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

    template <typename PathFollowingImpl, typename NonLinearOperator>
    concept PathFollowingStrategy
        = requires(PathFollowingImpl pft, NonLinearOperator nop, Ikarus::SubsidiaryArgs args) {
      { pft.evaluateSubsidiaryFunction(args) } -> std::same_as<void>;
      { pft.initialPrediction(nop, args) } -> std::same_as<void>;
      { pft.intermediatePrediction(nop, args) } -> std::same_as<void>;
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

    template <typename Op, typename... Args>
    concept IsFunctorWithArgs = requires(Op op, Args... args) {
      op(args...);
    };

    template <typename V>
    concept EigenVector = static_cast<bool>(V::IsVectorAtCompileTime);

#define MAKE_EIGEN_FIXED_VECTOR_CONCEPT(Size) \
  template <typename V>                       \
  concept EigenVector##Size                   \
      = static_cast<bool>(V::IsVectorAtCompileTime) and static_cast<bool>(V::SizeAtCompileTime == Size);

    MAKE_EIGEN_FIXED_VECTOR_CONCEPT(1)
    MAKE_EIGEN_FIXED_VECTOR_CONCEPT(2)
    MAKE_EIGEN_FIXED_VECTOR_CONCEPT(3)
    MAKE_EIGEN_FIXED_VECTOR_CONCEPT(4)
    MAKE_EIGEN_FIXED_VECTOR_CONCEPT(5)
    MAKE_EIGEN_FIXED_VECTOR_CONCEPT(6)

#define MAKE_EIGEN_FIXED_MATRIX_CONCEPT(Size1, Size2)                                                       \
  template <typename M>                                                                                     \
  concept EigenMatrix##Size1##Size2 = static_cast<bool>(std::remove_cvref_t<M>::RowsAtCompileTime == Size1) \
                                      and static_cast<bool>(std::remove_cvref_t<M>::ColsAtCompileTime == Size2);

    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 1)

    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 2)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 2)

    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 3)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 3)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 2)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 3)

    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 4)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 4)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 4)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 2)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 3)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 4)

    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 5)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 5)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 5)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 5)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 2)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 3)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 4)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 5)

    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(1, 6)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(2, 6)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(3, 6)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(4, 6)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(5, 6)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 2)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 3)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 4)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 5)
    MAKE_EIGEN_FIXED_MATRIX_CONCEPT(6, 6)

#define MAKE_EIGEN_FIXED_MATRIX_OR_VOIGT_CONCEPT(Size1, Size2) \
  template <typename M>                                        \
  concept EigenMatrixOrVoigtNotation##Size1 = EigenMatrix##Size1##Size1<M> or EigenVector##Size2<M>;

    MAKE_EIGEN_FIXED_MATRIX_OR_VOIGT_CONCEPT(1, 1)
    MAKE_EIGEN_FIXED_MATRIX_OR_VOIGT_CONCEPT(2, 3)
    MAKE_EIGEN_FIXED_MATRIX_OR_VOIGT_CONCEPT(3, 6)

    template <template <typename...> class MaterialToCheck, typename Material>
    consteval bool isMaterial() {
      if constexpr (Std::isSpecialization<MaterialToCheck, Material>::value) return true;

      if constexpr (Std::IsSpecializationNonTypeAndTypes<VanishingStress, Material>::value) {
        if constexpr (Std::isSpecialization<MaterialToCheck, typename Material::Underlying>::value) {
          return true;
        } else {
          return false;
        }
      } else
        return false;
    }

    template <template <typename...> class MaterialToCheck, typename Material>
    concept IsMaterial = isMaterial<MaterialToCheck, Material>();

  }  // namespace Concepts
}  // namespace Ikarus
