// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vanishingstrain.hh
 * \brief Defines the PlaneStrain material model and related functions.
 * \ingroup  materials
 */

// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/strainconversions.hh>
#include <ikarus/finiteelements/mechanics/materials/vanishingstress.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/nonlinearoperator.hh>

namespace Ikarus {

/**
 * \brief PlaneStrain material model that enforces strain components to be zero.
 * \ingroup materials
 * \tparam MI The underlying material model.
 */
template <auto strainIndexPair, typename MI>
struct VanishingStrain : public Material<VanishingStrain<strainIndexPair, MI>>
{
  /**
   * \brief Constructor for PlaneStrain.
   * \param mat The underlying material model.
   */
  explicit VanishingStrain(MI mat)
      : matImpl_{mat} {}

  using Underlying = MI;                              ///< The underlying material type.
  using ScalarType = typename Underlying::ScalarType; ///< Scalar type.

  [[nodiscard]] constexpr static std::string nameImpl() noexcept {
    auto matName = MI::name() + "_PlaneStrain";
    return matName;
  }
  static constexpr auto fixedPairs = strainIndexPair; ///< Array of fixed stress components.

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
   * \param E The Green-Lagrangian strain.
   * \return TangentModuli The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto Esol = reduceStrain(E);
    auto C          = matImpl_.template tangentModuli<Underlying::strainTag, true>(Esol);
    if constexpr (voigt)
      return staticCondensation(C, fixedVoigtIndices);
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

  /**
   * \brief Converts the input strain matrix to the appropriate form for stress reduction.
   * \tparam Derived The derived type of the input matrix.
   * \param E The input strain matrix.
   * \return decltype(auto) The converted strain matrix.
   */
  template <typename Derived>
  decltype(auto) maybeFromVoigt(const Eigen::MatrixBase<Derived>& E) const {
    if constexpr (Concepts::EigenVector<Derived>) { // receiving vector means Voigt notation
      return fromVoigt(E.derived(), true);
    } else
      return E.derived();
  }

  template <typename Derived>
  requires(strainTag != StrainTags::linear)
  auto reduceStrain(const Eigen::MatrixBase<Derived>& Eraw) const {
    decltype(auto) E                     = maybeFromVoigt(Eraw);
    std::remove_cvref_t<decltype(E)> Egl = transformStrain<strainTag, StrainTags::greenLagrangian>(E);

    for (auto [i, j] : fixedPairs) {
      Egl(i, j) = 0;
      Egl(j, i) = 0;
    }

    return transformStrain<StrainTags::greenLagrangian, strainTag>(Egl).derived();
  }

  template <typename Derived>
  requires(strainTag == StrainTags::linear)
  auto reduceStrain(const Eigen::MatrixBase<Derived>& Eraw) const {
    Eigen::Matrix3<ScalarType> E = maybeFromVoigt(Eraw);

    for (auto [i, j] : fixedPairs) {
      E(i, j) = 0;
      E(j, i) = 0;
    }

    return E;
  }
};

/**
 * \brief Factory function to create a VanishingStress material with specified stress indices.
 * \tparam stressIndexPair The array of StressIndexPair representing fixed stress components.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param p_tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material.
 */
template <Impl::StressIndexPair... stressIndexPair, typename MaterialImpl>
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
  return makeVanishingStrain<Impl::StressIndexPair{2, 1}, Impl::StressIndexPair{2, 0}, Impl::StressIndexPair{2, 2}>(
      mat);
}
} // namespace Ikarus