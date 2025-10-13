// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file decomposehyperelastic.hh
 * \brief Definition of helper struct and function to decompose a hyperelastic material model.
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus {
namespace Materials {

  /**
   * \brief A helper struct to access the return types of the function decomposeHyperelasticAndGetMaterialParameters.
   */
  template <typename MAT>
  struct DecomposedMaterialTypes
  {
    using TupleType = decltype(decomposeHyperelasticAndGetMaterialParameters(std::declval<MAT>()));
    using DevType   = std::tuple_element_t<0, TupleType>; ///< Type of the pure deviatoric material model.
    using VolType   = std::tuple_element_t<1, TupleType>; ///< Type of the pure volumetric material model.
    using DevParams = std::tuple_element_t<2, TupleType>; ///< Type of the parameters of the deviatoric material model.
    using VolParams = std::tuple_element_t<3, TupleType>; ///< Type of the parameters of the volumetric material model.
  };

  /**
   * \brief A helper function to decompose a hyperelastic material model and get the underlying deviatoric material
   * model, volumetric material model and their material parameters as a tuple. Here, the material models returned still
   * fulfill Concepts::Material and have also been reduced again if the underlying material model is reduced (vanishing
   * stress or vanishing strain).
   *
   * \remark It is to be noted that the volumetric material model returned comes with a material parameter equal to
   * unity. Hence, the volumetric material parameter of the underlying model is also returned separately. This is done
   * to obtain a simple constructor for the DisplacementPressure element.
   *
   * \tparam MAT Type of the material.
   * \param mat The underlying material model.
   * \return A tuple containing the deviatoric material model, volumetric material model, deviatoric material parameters
   * and volumetric material parameters.
   */
  template <Concepts::Material MAT>
  requires(MAT::isHyperelastic)
  auto decomposeHyperelasticAndGetMaterialParameters(const MAT& mat) {
    static_assert(
        Materials::UnderlyingMaterial_t<MAT>::hasVolumetricPart,
        "decomposeHyperelastic is only implemented for the case where both deviatoric and volumetric part exists.");
    const auto& umat                  = Materials::underlyingMaterial(mat);
    const auto dev                    = umat.deviatoricFunction();
    const auto vol                    = umat.volumetricFunction().volumetricFunction();
    const auto [devParams, volParams] = umat.materialParametersImpl();
    const auto hdev                   = Materials::Hyperelastic(dev);
    const auto hvol                   = Materials::makePureVolumetric(vol, 1.0); // material parameter set to unity.
    if constexpr (not MAT::isReduced) {
      return std::make_tuple(hdev, hvol, devParams, volParams);
    } else {
      constexpr auto fixedPairs = MAT::fixedPairs;
      if constexpr (MAT::isStrainVanished) {
        const auto reducedDev = Materials::VanishingStrain<fixedPairs, std::remove_cvref_t<decltype(hdev)>>(hdev);
        const auto reducedVol = Materials::VanishingStrain<fixedPairs, std::remove_cvref_t<decltype(hvol)>>(hvol);
        return std::make_tuple(reducedDev, reducedVol, devParams, volParams);
      } else {
        const auto reducedDev = Materials::VanishingStress<fixedPairs, std::remove_cvref_t<decltype(hdev)>>(hdev);
        const auto reducedVol = Materials::VanishingStress<fixedPairs, std::remove_cvref_t<decltype(hvol)>>(hvol);
        return std::make_tuple(reducedDev, reducedVol, devParams, volParams);
      }
    }
  }
} // namespace Materials
} // namespace Ikarus
