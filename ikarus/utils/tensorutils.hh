// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include <numeric>
#include <ranges>
#include <unsupported/Eigen/CXX11/Tensor>

#include <dune/common/promotiontraits.hh>

#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/math.hh>
namespace Ikarus {

  template <typename Derived, typename T, auto rank>
  Eigen::Tensor<typename Derived::Scalar, rank> TensorCast(const Eigen::EigenBase<Derived>& matrix,
                                                           const std::array<T, rank>& dims) {
    return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().eval().data(),
                                                                                       dims);
  }

  auto dyadic(const auto& A_ij, const auto& B_kl) {
    Eigen::array<Eigen::IndexPair<long>, 0> empty_index_list = {};
    return A_ij.contract(B_kl, empty_index_list).eval();
  }

  template <typename ScalarType = double, int dim = 3>
  auto symmetricIdentityFourthOrder() {
    Eigen::TensorFixedSize<double, Eigen::Sizes<dim, dim, dim, dim>> idTensor;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l)
            idTensor(i, j, k, l) = 0.5 * ((i == k) * (j == l) + (i == l) * (j == k));
    return idTensor;
  }

  template <typename ScalarType = double, int dim = 3>
  auto symmetricFourthOrder(const auto& A, const auto& B) {
    Eigen::TensorFixedSize<double, Eigen::Sizes<dim, dim, dim, dim>> idTensor;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l)
            idTensor(i, j, k, l) = 0.5 * (A(i, k) * B(j, l) + A(i, l) * B(j, k));
    return idTensor;
  }

  template <typename ScalarType = double, int dim = 3>
  auto identityFourthOrder() {
    Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>> idTensor;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l)
            idTensor(i, j, k, l) = (i == j) * (k == l);
    return idTensor;
  }

  template <typename AType, typename BType>
  auto fourthOrderIKJL(const Eigen::MatrixBase<AType>& A, const Eigen::MatrixBase<BType>& B) {
    static_assert(AType::RowsAtCompileTime == BType::RowsAtCompileTime);
    static_assert(AType::ColsAtCompileTime == BType::ColsAtCompileTime);
    using ScalarResultType =
        typename Dune::PromotionTraits<typename AType::Scalar, typename BType::Scalar>::PromotedType;
    constexpr int dim = AType::RowsAtCompileTime;
    Eigen::TensorFixedSize<ScalarResultType, Eigen::Sizes<dim, dim, dim, dim>> res;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l)
            res(i, j, k, l) = A(i, k) * B(j, l);
    return res;
  }

  template <typename ScalarType, long int dim>
  auto symTwoSlots(const Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>>& t,
                   const std::array<size_t, 2>& slots) {
    std::array<size_t, 4> shuffleSlot;
    std::iota(shuffleSlot.begin(), shuffleSlot.end(), 0);     // create 0,1,2,3 array
    std::swap(shuffleSlot[slots[0]], shuffleSlot[slots[1]]);  // swap  the given slots
    return (0.5 * (t + t.shuffle(shuffleSlot))).eval();
  }

  constexpr Eigen::Index toVoigt(Eigen::Index i, Eigen::Index j) noexcept {
    if (i == j)  // _00 -> 0, _11 -> 1,  _22 -> 2
      return i;
    else if ((i == 1 and j == 2) or (i == 2 and j == 1))  // _12 and _21 --> 3
      return 3;
    else if ((i == 0 and j == 2) or (i == 2 and j == 0))  // _02 and _20 --> 4
      return 4;
    else if ((i == 0 and j == 1) or (i == 1 and j == 0))  // _01 and _10 --> 5
      return 5;
    else {
      assert(i < 3 and j < 3 && "For Voigt notation the indices need to be 0,1 or 2.");
      __builtin_unreachable();
    }
  }

  template <typename ScalarType = double>
  Eigen::Matrix<ScalarType, 6, 6> toVoigt(const Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>& ft) {
    Eigen::Matrix<ScalarType, 6, 6> mat;
    for (Eigen::Index i = 0; i < 3; ++i)
      for (Eigen::Index j = 0; j < 3; ++j)
        for (Eigen::Index k = 0; k < 3; ++k)
          for (Eigen::Index l = 0; l < 3; ++l)
            mat(toVoigt(i, j), toVoigt(k, l)) = ft(i, j, k, l);  // TODO exploit symmetry
    return mat;
  }

  template <typename ST, int size, int Options, int MaxRows>
  requires(size > 0 and size <= 3) auto toVoigt(const Eigen::Matrix<ST, size, size, Options, MaxRows, MaxRows>& E,
                                                bool isStrain = true) {
    Eigen::Vector<ST, (size * (size + 1)) / 2> EVoigt;
    EVoigt.template segment<size>(0) = E.diagonal();

    const ST possibleStrainFactor = isStrain ? 2.0 : 1.0;
    if constexpr (size == 2)
      EVoigt(2) = E(0, 1) * possibleStrainFactor;
    else if constexpr (size == 3) {
      EVoigt(size)     = E(1, 2) * possibleStrainFactor;
      EVoigt(size + 1) = E(0, 2) * possibleStrainFactor;
      EVoigt(size + 2) = E(0, 1) * possibleStrainFactor;
    }
    return EVoigt;
  }

  constexpr std::array<size_t, 2> fromVoigt(size_t i) {
    if (i < 3)  // _00 -> 0, _11 -> 1,  _22 -> 2
      return {i, i};
    else if (i == 3)
      return {1, 2};
    else if (i == 4)
      return {0, 2};
    else if (i == 5)
      return {0, 1};
    else {
      assert(i < 6 && "For Voigt notation the indices need to be 0 and 5.");
      __builtin_unreachable();
    }
  }

  template <typename ST, int size>
  requires(size == 1 or size == 3 or size == 6) auto fromVoigt(const Eigen::Vector<ST, size>& EVoigt,
                                                               bool isStrain = true) {
    constexpr int matrixSize = (-1 + Ikarus::ct_sqrt(1 + 8 * size)) / 2;
    Eigen::Matrix<ST, matrixSize, matrixSize> E;
    E.diagonal() = EVoigt.template segment<3>(0);

    const ST possibleStrainFactor = isStrain ? 0.5 : 1.0;
    if constexpr (matrixSize == 2)
      E(0, 1) = E(1, 0) = EVoigt(2) * possibleStrainFactor;
    else if constexpr (matrixSize == 3) {
      E(2, 1) = E(1, 2) = EVoigt(matrixSize) * possibleStrainFactor;
      E(2, 0) = E(0, 2) = EVoigt(matrixSize + 1) * possibleStrainFactor;
      E(1, 0) = E(0, 1) = EVoigt(matrixSize + 2) * possibleStrainFactor;
    }
    return E;
  }

  template <typename ScalarType>
  auto fromVoigt(const Eigen::Matrix<ScalarType, 6, 6>& CVoigt) {
    Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> C;
    //    size_t iR=0,jR=0;
    for (size_t i = 0; i < 6; ++i) {
      for (size_t j = 0; j < 6; ++j) {
        auto firstIndices                                                       = fromVoigt(i);
        auto secondIndices                                                      = fromVoigt(j);
        C(firstIndices[0], firstIndices[1], secondIndices[0], secondIndices[1]) = CVoigt(i, j);
      }  // TODO exploit symmetry
    }
    return C;
  }

  template <typename Derived, size_t sizeOfCondensedIndices>
  auto staticCondensation(const Eigen::MatrixBase<Derived>& E,
                          const std::array<size_t, sizeOfCondensedIndices>& indices) {
    constexpr size_t colsFull = std::remove_cvref_t<Derived>::ColsAtCompileTime;
    constexpr size_t rowsFull = std::remove_cvref_t<Derived>::RowsAtCompileTime;
    static_assert(colsFull == rowsFull, "staticCondensation only allowed for square matrices");
    std::array<size_t, rowsFull - sizeOfCondensedIndices> remainingIndices{};
    std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(colsFull)), indices, remainingIndices.begin());

    auto K11 = E(remainingIndices, remainingIndices);
    auto K12 = E(indices, remainingIndices);
    auto K22 = E(indices, indices);

    return (K11 - K12.transpose() * K22.inverse() * K12).eval();
  }

  template <typename Derived, size_t sizeOfRemovedCols>
  auto removeCol(const Eigen::MatrixBase<Derived>& E, const std::array<size_t, sizeOfRemovedCols>& indices) {
    constexpr size_t colsFull = std::remove_cvref_t<Derived>::ColsAtCompileTime;
    constexpr size_t rowsFull = std::remove_cvref_t<Derived>::RowsAtCompileTime;
    static_assert(colsFull == 1);

    std::array<size_t, rowsFull - sizeOfRemovedCols> remainingIndices{};
    std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(rowsFull)), indices, remainingIndices.begin());

    return (E(remainingIndices)).eval();
  }

  template <typename ST, typename MaterialImpl>
  auto toVoigtAndMaybeReduce(const Eigen::Matrix<ST, 3, 3>& E, const MaterialImpl&, bool isStrain = true) {
    if constexpr (!MaterialImpl::isReduced)
      return toVoigt(E, isStrain);
    else {
      auto ev = toVoigt(E, isStrain);
      static_assert(decltype(ev)::RowsAtCompileTime == 6 and decltype(ev)::ColsAtCompileTime == 1);

      auto evRed = removeCol(ev, MaterialImpl::fixedVoigtIndices);
      static_assert(decltype(evRed)::RowsAtCompileTime == 6 - MaterialImpl::fixedVoigtIndices.size()
                    and decltype(evRed)::ColsAtCompileTime == 1);
      return evRed;
    }
  }

  template <typename Material, typename Derived>
  decltype(auto) enlargeIfReduced(const Eigen::MatrixBase<Derived>& E) {
    if constexpr (!Concepts::EigenVector6<Derived> and Concepts::EigenVector<Derived>) {
      static_assert(Material::isReduced, "You should only end up here, if your material is reduced");

      auto freeindices = Material::MaterialImpl::freeVoigtIndices;
      auto p_E         = Eigen::Vector<typename Material::MaterialImpl::ScalarType, 6>::Zero().eval();
      for (int ri = 0; auto i : freeindices)
        p_E(i) = E(ri++);
      return p_E;

    } else if constexpr (Concepts::EigenMatrix66<
                             Derived> or Concepts::EigenMatrix33<Derived> or Concepts::EigenVector6<Derived>) {
      return E.derived();
    } else {
      static_assert(Material::isReduced, "You should only end up here, if your material is reduced");

      auto freeindices = Material::MaterialImpl::freeVoigtIndices;
      auto p_E         = Eigen::Matrix<typename Material::MaterialImpl::ScalarType, 6, 6>::Zero().eval();
      for (int ri = 0; auto i : freeindices) {
        for (int rj = 0; auto j : freeindices)
          p_E(i, j) = E(ri, rj++);
        ++ri;
      }
      return p_E;
    }
  }

}  // namespace Ikarus
