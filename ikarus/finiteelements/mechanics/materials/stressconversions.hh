// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file stressconversions.hh
 * \brief Implementation of transformations for different stress measures
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
 * \brief Create PK2 stress matrix based on the input.
 *
 * This function creates PK2 stress matrix based on the input stress matrix and deformation gradient.
 * What to do is decided by the provided stress tag
 * \tparam tag Tag of the input stress measure.
 * \tparam Derived Type of the Eigen matrices.
 * \param sMB Input stress matrix.
 * \param F The deformation gradient.
 * \return The PK2 stress matrix.
 */
template <StressTags tag, typename Derived>
Derived createPK2Stress(const Eigen::MatrixBase<Derived>& sMB, const Eigen::MatrixBase<Derived>& F) {
  const auto& S = sMB.derived();

  static_assert(Concepts::EigenMatrix33<Derived> or Concepts::EigenMatrix22<Derived>);
  if constexpr (tag == StressTags::Cauchy) {
    auto J    = F.determinant();
    auto invF = F.inverse().eval();
    return (J * invF * S * invF.transpose()).eval();
  } else if constexpr (tag == StressTags::Kirchhoff) {
    auto invF = F.inverse().eval();
    return (invF * S * invF.transpose()).eval();
  } else if constexpr (tag == StressTags::PK1)
    return (F.inverse() * S).eval();
  else if constexpr (tag == StressTags::PK2)
    return S;
}

/**
 * \brief Create PK1 stress matrix based on the input.
 *
 * This function creates PK1 stress matrix based on the input stress matrix and deformation gradient.
 * What to do is decided by the provided stress tag
 * \tparam tag Tag of the input stress measure.
 * \tparam Derived Type of the Eigen matrices.
 * \param sMB Input stress matrix.
 * \param F The deformation gradient.
 * \return The PK1 stress matrix.
 */
template <StressTags tag, typename Derived>
Derived createPK1Stress(const Eigen::MatrixBase<Derived>& sMB, const Eigen::MatrixBase<Derived>& F) {
  const auto& S = sMB.derived();

  static_assert(Concepts::EigenMatrix33<Derived> or Concepts::EigenMatrix22<Derived>);
  if constexpr (tag == StressTags::Cauchy)
    return (F.determinant() * S * F.inverse().transpose()).eval();
  else if constexpr (tag == StressTags::Kirchhoff)
    return (S * F.inverse().transpose()).eval();
  else if constexpr (tag == StressTags::PK1)
    return S;
  else if constexpr (tag == StressTags::PK2)
    return (F * S).eval();
}

/**
 * \brief Create Kirchhoff stress matrix based on the input.
 *
 * This function creates Kirchhoff stress matrix based on the input stress matrix and deformation gradient.
 * What to do is decided by the provided stress tag
 * \tparam tag Tag of the input stress measure.
 * \tparam Derived Type of the Eigen matrices.
 * \param sMB Input stress matrix.
 * \param F The deformation gradient.
 * \return The Kirchhoff stress matrix.
 */
template <StressTags tag, typename Derived>
Derived createKirchhoffStress(const Eigen::MatrixBase<Derived>& sMB, const Eigen::MatrixBase<Derived>& F) {
  const auto& S = sMB.derived();

  static_assert(Concepts::EigenMatrix33<Derived> or Concepts::EigenMatrix22<Derived>);
  if constexpr (tag == StressTags::Cauchy)
    return (F.determinant() * S).eval();
  else if constexpr (tag == StressTags::Kirchhoff)
    return S;
  else if constexpr (tag == StressTags::PK1)
    return (S * F.transpose()).eval();
  else if constexpr (tag == StressTags::PK2)
    return (F * S * F.transpose()).eval();
}

/**
 * \brief Create Cauchy stress matrix based on the input.
 *
 * This function creates Cauchy stress matrix based on the input stress matrix and deformation gradient.
 * What to do is decided by the provided stress tag
 * \tparam tag Tag of the input stress measure.
 * \tparam Derived Type of the Eigen matrices.
 * \param sMB Input stress matrix.
 * \param F The deformation gradient.
 * \return The Cauchy stress matrix.
 */
template <StressTags tag, typename Derived>
Derived createCauchyStress(const Eigen::MatrixBase<Derived>& sMB, const Eigen::MatrixBase<Derived>& F) {
  const auto& S   = sMB.derived();
  const auto invJ = 1.0 / F.determinant();

  static_assert(Concepts::EigenMatrix33<Derived> or Concepts::EigenMatrix22<Derived>);
  if constexpr (tag == StressTags::Cauchy)
    return S;
  else if constexpr (tag == StressTags::Kirchhoff)
    return (invJ * S).eval();
  else if constexpr (tag == StressTags::PK1)
    return (invJ * S * F.transpose()).eval();
  else if constexpr (tag == StressTags::PK2)
    return (invJ * F * S * F.transpose()).eval();
}

/**
 * \brief Transform stress measures from one type to another.
 *
 * This function transforms a stress component matrix from one type to another, based on the provided stress tags and
 * the deformation gradient
 * \tparam from Tag of the source stress measure.
 * \tparam to Tag of the target stress measure.
 * \tparam DerivedS Type of the stress matrix.
 * \param sRaw Eigen matrix representing the input stress (can be in Voigt notation).
 * \param F Eigen matrix representing the deformation gradient (has to be in matrix notation).
 * \return The transformed stress matrix.
 */
template <StressTags from, StressTags to, typename DerivedS, typename DerivedF>
auto transformStress(const Eigen::MatrixBase<DerivedS>& sRaw, const Eigen::MatrixBase<DerivedF>& F) {
  static_assert((from == to) or (from != StressTags::linear and to != StressTags::linear),
                "No useful transformation available for linear stresses.");
  static_assert(Concepts::EigenMatrix33<DerivedF> or Concepts::EigenMatrix22<DerivedF>);

  const auto S = Impl::maybeFromVoigt(sRaw.derived(), false);

  if constexpr (from == to)
    return S;
  else if constexpr (to == StressTags::PK2)
    return createPK2Stress<from>(S, F);
  else if constexpr (to == StressTags::PK1)
    return createPK1Stress<from>(S, F);
  else if constexpr (to == StressTags::Kirchhoff)
    return createKirchhoffStress<from>(S, F);
  else if constexpr (to == StressTags::Cauchy)
    return createCauchyStress<from>(S, F);
  else
    static_assert(to == StressTags::PK2 or to == StressTags::PK1 or to == StressTags::Kirchhoff or
                  to == StressTags::Cauchy);
}
} // namespace Ikarus
