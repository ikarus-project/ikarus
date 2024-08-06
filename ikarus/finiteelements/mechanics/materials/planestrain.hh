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
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/nonlinearoperator.hh>

namespace Ikarus {

/**
 * \brief PlaneStrain material model that enforces strain components to be zero.
 * \ingroup materials
 * \tparam MI The underlying material model.
 */
template <typename MI>
struct PlaneStrain : public Material<PlaneStrain<MI>>
{
  /**
   * \brief Constructor for PlaneStrain.
   * \param mat The underlying material model.
   */
  explicit PlaneStrain(MI mat)
      : matImpl_{mat} {}

  using Underlying = MI;                              ///< The underlying material type.
  using ScalarType = typename Underlying::ScalarType; ///< Scalar type.

  [[nodiscard]] constexpr std::string nameImpl() const noexcept {
    auto matName = matImpl_.name() + "_PlaneStrain";
    return matName;
  }

  static constexpr auto strainTag          = Underlying::strainTag;          ///< Strain tag.
  static constexpr auto stressTag          = Underlying::stressTag;          ///< Stress tag.
  static constexpr auto tangentModuliTag   = Underlying::tangentModuliTag;   ///< Tangent moduli tag.
  static constexpr bool energyAcceptsVoigt = Underlying::energyAcceptsVoigt; ///< Energy accepts Voigt notation.
  static constexpr bool stressToVoigt      = true;                           ///< Stress to Voigt notation.
  static constexpr bool stressAcceptsVoigt = true;                           ///< Stress accepts Voigt notation.
  static constexpr bool moduliToVoigt      = true;                           ///< Moduli to Voigt notation.
  static constexpr bool moduliAcceptsVoigt = true;                           ///< Moduli accepts Voigt notation.

  static constexpr auto freeStrains       = 3;                              ///< Number of free strains.
  static constexpr auto freeVoigtIndices  = std::array<size_t, 3>{0, 1, 5}; ///< Free Voigt indices.
  static constexpr auto fixedVoigtIndices = std::array<size_t, 3>{2, 3, 4}; ///< Fixed Voigt indices.

  /**
   * \brief Computes the stored energy for the PlaneStrain material.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
    return matImpl_.storedEnergyImpl(E);
  }

  /**
   * \brief Computes the straines for the PlaneStrain material.
   * \tparam voigt A boolean indicating whether to return straines in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return StressMatrix The straines.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& Eraw) const {
    auto E         = maybeFromVoigt(Eraw);
    auto stresses  = matImpl_.template stresses<Underlying::strainTag, true>(E);
    auto stressRed = stresses({0, 1, 5}).eval();
    if constexpr (voigt)
      return stressRed;
    else
      return fromVoigt(stressRed, false);
  }

  /**
   * \brief Computes the tangent moduli for the PlaneStrain material.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return TangentModuli The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& Eraw) const {
    auto E                              = maybeFromVoigt(Eraw);
    auto C                              = matImpl_.template tangentModuli<Underlying::strainTag, true>(E);
    Eigen::Matrix<ScalarType, 3, 3> C33 = C({0, 1, 5}, {0, 1, 5}).eval();
    if constexpr (voigt)
      return C33;
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
};

/**
 * \brief Factory function to create a PlaneStrain material.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \return PlaneStrain The created PlaneStrain material.
 */
template <typename MaterialImpl>
auto makePlaneStrain(MaterialImpl mat) {
  return PlaneStrain<MaterialImpl>(mat);
}

/**
 * \brief Factory function to create a PlaneStrain material for plane strain conditions.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \return PlaneStrain The created PlaneStrain material for plane strain case.
 */
template <typename MaterialImpl>
auto planeStrain(const MaterialImpl& mat) {
  return makePlaneStrain(mat);
}
} // namespace Ikarus
