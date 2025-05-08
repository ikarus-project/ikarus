// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file strainconversions.hh
 * \brief Implementation of transformations for different strain measures
 *
 * \ingroup materials
 */

#pragma once
#include "tags.hh"

#include <unsupported/Eigen/MatrixFunctions>

#include <Eigen/Core>

#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

/**
 * \brief Create Green-Lagrangian strain based on the input.
 *
 * This function creates Green-Lagrangian strains based on the input strain matrix.
 * What to do is decided by the provided strain tag
 * \tparam tag Tag of the input strain measure.
 * \tparam Derived Type of the Eigen matrix.
 * \param eMB Eigen matrix representing the input strain.
 * \return The Green-Lagrangian strains matrix.
 */
template <StrainTags tag, typename Derived>
Derived createGreenLagrangianStrains(const Eigen::MatrixBase<Derived>& eMB) {
  const auto& e = eMB.derived();

  static_assert(Concepts::EigenMatrix33<Derived> or Concepts::EigenMatrix22<Derived>);
  if constexpr (tag == StrainTags::greenLagrangian)
    return e;
  else if constexpr (tag == StrainTags::deformationGradient)
    return 0.5 * (e.transpose() * e - Derived::Identity());
  else if constexpr (tag == StrainTags::displacementGradient)
    return 0.5 * (e + e.transpose() + e.transpose() * e);
  else if constexpr (tag == StrainTags::rightCauchyGreenTensor)
    return 0.5 * (e - Derived::Identity());
}

/**
 * \brief Create the deformation gradient based on the input.
 *
 * This function creates deformation gradient based on the input strain matrix.
 * What to do is decided by the provided strain tag
 * \tparam tag Tag of the input strain measure.
 * \tparam Derived Type of the Eigen matrix.
 * \param eMB Eigen matrix representing the input strain.
 * \return The deformation gradient matrix.
 */
template <StrainTags tag, typename Derived>
Derived createDeformationGradient(const Eigen::MatrixBase<Derived>& eMB) {
  const auto& e = eMB.derived();

  static_assert(Concepts::EigenMatrix33<Derived> or Concepts::EigenMatrix22<Derived>);
  if constexpr (tag == StrainTags::greenLagrangian) {
    // E = 0.5 * (F ^ 2 - I);
    // 2*E = F ^ 2 - I;
    // 2*E+I = F ^ 2;
    // sqrt(2*E+I) = F;
    return (2 * e + Derived::Identity()).sqrt();
  } else if constexpr (tag == StrainTags::deformationGradient)
    return e;
  else if constexpr (tag == StrainTags::displacementGradient)
    return e + Derived::Identity();
  else if constexpr (tag == StrainTags::rightCauchyGreenTensor) {
    return e.sqrt(); // this looses information, since the rotation information from the original F is lost!
  }
}

/**
 * \brief Create right Cauchy-Green tensor based on the input.
 *
 * This function creates Right Cauchy-Green tensor based on the input strain matrix.
 * What to do is decided by the provided strain tag
 * \tparam tag Tag of the input strain measure.
 * \tparam Derived Type of the Eigen matrix.
 * \param eMB Eigen matrix representing the input strain.
 * \return The Right Cauchy-Green tensor matrix.
 */
template <StrainTags tag, typename Derived>
Derived createRightCauchyGreen(const Eigen::MatrixBase<Derived>& eMB) {
  const auto& e = eMB.derived();

  static_assert(Concepts::EigenMatrix33<Derived> or Concepts::EigenMatrix22<Derived>);
  if constexpr (tag == StrainTags::greenLagrangian) {
    // E = 0.5 * (C - I);
    // 2*E = C - I;
    // 2*E+I = C;
    return 2 * e + Derived::Identity();
  } else if constexpr (tag == StrainTags::deformationGradient)
    return e.transpose() * e;
  else if constexpr (tag == StrainTags::displacementGradient) {
    const auto F = e + Derived::Identity();
    return F.transpose() * F;
  } else if constexpr (tag == StrainTags::rightCauchyGreenTensor) {
    return e;
  }
}

/**
 * \brief Transform strain from one type to another.
 *
 * This function transforms one strain component matrix from one type to another, based on the provided strain tags
 * \tparam from Tag of the source strain measure.
 * \tparam to Tag of the target strain measure.
 * \tparam Derived Type of the Eigen matrix.
 * \param eRaw Eigen matrix representing the input strain (can be in Voigt notation).
 * \return The transformed strain matrix.
 */
template <StrainTags from, StrainTags to, typename Derived>
auto transformStrain(const Eigen::MatrixBase<Derived>& eRaw) {
  static_assert((from == to) or (from != StrainTags::linear and to != StrainTags::linear),
                "No useful transformation available for linear strains.");
  const auto e = Impl::maybeFromVoigt(eRaw.derived(), true);
  if constexpr (from == to)
    return e;
  else if constexpr (to == StrainTags::greenLagrangian)
    return createGreenLagrangianStrains<from>(e);
  else if constexpr (to == StrainTags::deformationGradient)
    return createDeformationGradient<from>(e);
  else if constexpr (to == StrainTags::rightCauchyGreenTensor) {
    return createRightCauchyGreen<from>(e);
  } else
    static_assert(to == StrainTags::greenLagrangian or to == StrainTags::deformationGradient or
                  to == StrainTags::rightCauchyGreenTensor);
}
} // namespace Ikarus
