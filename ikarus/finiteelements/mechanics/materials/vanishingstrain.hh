// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vanishingstrain.hh
 * \brief Defines the VanishingStrain material model and related functions.
 * \ingroup  materials
 */

// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "vanishinghelpers.hh"

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/strainconversions.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/nonlinearoperator.hh>

namespace Ikarus {

/**
 * \brief VanishingStrain material model that enforces strain components to be zero.
 * \ingroup materials
 * \tparam strainIndexPair An array of MatrixIndexPair representing fixed stress components.
 * \tparam MI The underlying material model.
 */
template <auto strainIndexPair, typename MI>
struct VanishingStrain : public Material<VanishingStrain<strainIndexPair, MI>>
{
  using Underlying         = MI;                              ///< The underlying material type.
  using ScalarType         = typename Underlying::ScalarType; ///< Scalar type.
  using MaterialParameters = typename Underlying::MaterialParameters;

  static constexpr auto fixedPairs        = strainIndexPair;                     ///< Array of fixed stress components.
  static constexpr auto freeVoigtIndices  = createfreeVoigtIndices(fixedPairs);  ///< Free Voigt indices.
  static constexpr auto fixedVoigtIndices = createFixedVoigtIndices(fixedPairs); ///< Fixed Voigt indices.
  static constexpr auto freeStrains       = freeVoigtIndices.size();             ///< Number of free strains.

  static constexpr auto strainTag          = Underlying::strainTag;          ///< Strain tag.
  static constexpr auto stressTag          = Underlying::stressTag;          ///< Stress tag.
  static constexpr auto tangentModuliTag   = Underlying::tangentModuliTag;   ///< Tangent moduli tag.
  static constexpr bool energyAcceptsVoigt = Underlying::energyAcceptsVoigt; ///< Energy accepts Voigt notation.
  static constexpr bool stressToVoigt      = true;                           ///< Stress to Voigt notation.
  static constexpr bool stressAcceptsVoigt = true;                           ///< Stress accepts Voigt notation.
  static constexpr bool moduliToVoigt      = true;                           ///< Moduli to Voigt notation.
  static constexpr bool moduliAcceptsVoigt = true;                           ///< Moduli accepts Voigt notation.
  static constexpr double derivativeFactor = 1;                              ///< Derivative factor.

  /**
   * \brief Constructor for VanishingStrain.
   * \param mat The underlying material model.
   */
  explicit VanishingStrain(MI mat)
      : matImpl_{mat} {}

  [[nodiscard]] constexpr static std::string nameImpl() noexcept {
    auto matName = MI::name() + "_VanishingStrain(";
    for (auto p : fixedPairs)
      matName += "(" + std::to_string(p.row) + std::to_string(p.col) + ")";
    matName += ")";
    return matName;
  }

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return matImpl_.materialParametersImpl(); }

  /**
   * \brief Computes the stored energy for the PlaneStrain material.
   * \tparam Derived The derived type of the input matrix.
   * \param Eraw The strain mesasure
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto Esol = reduceStrain(E);
    return matImpl_.storedEnergyImpl(Esol);
  }

  /**
   * \brief Computes the strains for the PlaneStrain material.
   * \tparam voigt A boolean indicating whether to return strains in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param Eraw The Green-Lagrangian strain.
   * \return StressMatrix The strains.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto Esol  = reduceStrain(E);
    auto stressesRed = matImpl_.template stresses<Underlying::strainTag, true>(Esol);

    if constexpr (voigt) {
      return removeCol(stressesRed, fixedVoigtIndices);
    } else {
      stressesRed(fixedVoigtIndices).setZero();
      return fromVoigt(stressesRed, false);
    }
  }

  /**
   * \brief Computes the tangent moduli for the PlaneStrain material.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The strain measure.
   * \return TangentModuli The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto Esol = reduceStrain(E);
    auto C          = matImpl_.template tangentModuli<Underlying::strainTag, true>(Esol);
    if constexpr (voigt)
      return reduceMatrix(C, fixedVoigtIndices);
    else
      return fromVoigt(C);
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam ScalarTypeOther The target scalar type.
   * \return PlaneStrain The rebound PlaneStrain material.
   */
  template <typename ScalarTypeOther>
  auto rebind() const {
    auto reboundMatImpl = matImpl_.template rebind<ScalarTypeOther>();
    return PlaneStrain<decltype(reboundMatImpl)>(reboundMatImpl);
  }

private:
  Underlying matImpl_; ///< The underlying material model.

  template <typename Derived>
  auto reduceStrain(const Eigen::MatrixBase<Derived>& Eraw) const {
    if constexpr (strainTag == StrainTags::linear or strainTag == StrainTags::greenLagrangian) {
      Eigen::Matrix3<ScalarType> E = Impl::maybeFromVoigt(Eraw);
      setStrainsToZero(E);
      return E;
    } else {
      decltype(auto) E               = Impl::maybeFromVoigt(Eraw);
      Eigen::Matrix3<ScalarType> Egl = transformStrain<strainTag, StrainTags::greenLagrangian>(E);
      setStrainsToZero(Egl);
      return transformStrain<StrainTags::greenLagrangian, strainTag>(Egl).derived();
    }
  }

  inline void setStrainsToZero(auto& E) const {
    for (auto [i, j] : fixedPairs) {
      E(i, j) = 0;
      E(j, i) = 0;
    }
  }
};

/**
 * \brief Factory function to create a PlaneStrain material with specified strain indices.
 * \tparam matrixIndexPair The array of MatrixIndexPair representing fixed strain components.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param p_tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material.
 */
template <Impl::MatrixIndexPair... stressIndexPair, typename MaterialImpl>
auto makeVanishingStrain(MaterialImpl mat) {
  return VanishingStrain<std::to_array({stressIndexPair...}), MaterialImpl>(mat);
}

/**
 * \brief Factory function to create a PlaneStrain material for plane strain conditions.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \return PlaneStrain The created PlaneStrain material for plane strain case.
 */
template <typename MaterialImpl>
auto planeStrain(const MaterialImpl& mat) {
  return makeVanishingStrain<Impl::MatrixIndexPair{2, 1}, Impl::MatrixIndexPair{2, 0}, Impl::MatrixIndexPair{2, 2}>(
      mat);
}
} // namespace Ikarus