// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materialconversions.hh
 * \brief Implementation of transformations for different material tensor measures
 *
 * \ingroup materials
 */

#pragma once
#include "tags.hh"

#include <Eigen/Core>

#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

/**
 * \brief Create two-point material tensor based on the input.
 *
 * \details This function creates two-point material tensor based on the input material tensor, PK2 stress tensor and
 * deformation gradient. What to do is decided by the provided tangent moduli tag. Voigt notation is generally not
 * applicable for the two-point material tensor. The transformation from spatial to two-point tensor is not implemented.
 * All input quantities must be in 3D notations, i.e., no Voigt notation or reduction due to plane stress or plane
 * strain case is considered.
 * \tparam tag Tag of the input measure.
 * \tparam Derived Type of the Eigen matrices.
 * \tparam ST Scalar type of the entries of the input tangent moduli.
 * \param C Input tangent moduli tensor (in tensor notation)
 * \param sPK2 Input PK2 stress matrix (in matrix notation)
 * \param F The deformation gradient (in matrix notation)
 * \return The two-point material tensor (in tensor notation)
 */
template <TangentModuliTags tag, typename Derived, typename ST>
Eigen::TensorFixedSize<ST, Eigen::Sizes<3, 3, 3, 3>> createTwoPointMaterialTensor(
    const Eigen::TensorFixedSize<ST, Eigen::Sizes<3, 3, 3, 3>>& C, const Eigen::MatrixBase<Derived>& sPK2,
    const Eigen::MatrixBase<Derived>& F) {
  static_assert(Concepts::EigenMatrix33<Derived>);
  if constexpr (tag == TangentModuliTags::Material) {
    const auto& S                    = sPK2.derived();
    const Eigen::Matrix<ST, 3, 3> Id = Eigen::Matrix<ST, 3, 3>::Identity();
    constexpr int dim                = 3;
    Eigen::TensorFixedSize<ST, Eigen::Sizes<dim, dim, dim, dim>> A;
    A.setZero();
    for (const auto i : Dune::range(dim))
      for (const auto J : Dune::range(dim))
        for (const auto k : Dune::range(dim))
          for (const auto L : Dune::range(dim))
            for (const auto I : Dune::range(dim))
              for (const auto K : Dune::range(dim)) {
                A(i, J, k, L) += C(I, J, K, L) * F(i, I) * F(k, K) + Id(i, k) * S(J, L);
              }
    return A;
  } else if constexpr (tag == TangentModuliTags::TwoPoint) {
    return C;
  } else if constexpr (tag == TangentModuliTags::Spatial) {
    static_assert(Dune::AlwaysFalse<ST>::value,
                  "Transformation from spatial to two-point material tensor is not implemented.");
  }
}

/**
 * \brief Transform tangent moduli measures from one type to another.
 *
 * \details This function transforms a tangent moduli tensor from one type to another, based on the provided tangent
 * moduli tags, PK2 stress tensor and the deformation gradient
 * \tparam from Tag of the source tangent moduli measure.
 * \tparam to Tag of the target tangent moduli measure.
 * \tparam DerivedS Type of the stress matrix.
 * \tparam DerivedF Type of the deformation gradient matrix.
 * \tparam ST Scalar type of the entries of the input tangent moduli.
 * \param sRaw Eigen matrix representing the input PK2 stress (can be in Voigt notation).
 * \param F Eigen matrix representing the deformation gradient (can only be a 3D matrix).
 * \return The transformed tangent moduli in tensor notation.
 */
template <TangentModuliTags from, TangentModuliTags to, typename DerivedS, typename DerivedF, typename ST>
auto transformTangentModuli(const Eigen::TensorFixedSize<ST, Eigen::Sizes<3, 3, 3, 3>>& C,
                            const Eigen::MatrixBase<DerivedS>& sRaw, const Eigen::MatrixBase<DerivedF>& F) {
  static_assert(Concepts::EigenMatrix33<DerivedF>);

  const auto S = Impl::maybeFromVoigt(sRaw.derived(), false);

  if constexpr (from == to)
    return C;
  else if constexpr (to == TangentModuliTags::TwoPoint)
    return createTwoPointMaterialTensor<from>(C, S, F);
  else
    static_assert(Dune::AlwaysFalse<ST>::value, "Transformation only to two-point material tensor is implemented.");
}
} // namespace Ikarus
