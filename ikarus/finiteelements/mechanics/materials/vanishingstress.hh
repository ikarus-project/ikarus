// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vanishingstress.hh
 * \brief Defines the VanishingStress material model and related functions.
 * \ingroup  materials
 */

// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/nonlinearoperator.hh>

namespace Ikarus {

namespace Impl {

  /**
   * \brief Represents a pair of stress matrix indices (row and column).
   */
  struct StressIndexPair
  {
    Eigen::Index row; ///< Row index.
    Eigen::Index col; ///< Column index.
  };

  /**
   * \brief Helper function to create an array of free Voigt indices.
   * \tparam size The size of the fixed pairs array.
   * \param fixed An array of StressIndexPair representing fixed indices.
   * \return std::array<size_t, 6 - size> The array of free Voigt indices.
   */
  template <size_t size>
  consteval auto createfreeVoigtIndices(const std::array<StressIndexPair, size>& fixed) {
    std::array<size_t, 6 - size> res{};
    std::array<size_t, size> voigtFixedIndices;
    std::ranges::transform(fixed, voigtFixedIndices.begin(), [](auto pair) { return toVoigt(pair.row, pair.col); });
    std::ranges::sort(voigtFixedIndices);
    std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(6)), voigtFixedIndices, res.begin());
    std::ranges::sort(res);
    return res;
  }

  /**
   * \brief Helper function to create an array of fixed Voigt indices.
   * \tparam size The size of the fixed pairs array.
   * \param fixed An array of StressIndexPair representing fixed indices.
   * \return std::array<size_t, size> The array of fixed Voigt indices.
   */
  template <size_t size>
  consteval auto createFixedVoigtIndices(const std::array<StressIndexPair, size>& fixed) {
    std::array<size_t, size> fixedIndices;
    std::ranges::transform(fixed, fixedIndices.begin(), [](auto pair) { return toVoigt(pair.row, pair.col); });
    std::ranges::sort(fixedIndices);
    return fixedIndices;
  }

  /**
   * \brief Helper function to count the number of diagonal indices in the fixed pairs array.
   * \tparam size The size of the fixed pairs array.
   * \param fixed An array of StressIndexPair representing fixed indices.
   * \return constexpr size_t The number of diagonal indices.
   */
  template <size_t size>
  constexpr size_t countDiagonalIndices(const std::array<StressIndexPair, size>& fixed) {
    size_t count = 0;
    for (auto v : fixed) {
      if (v.col == v.row)
        ++count;
    }
    return count;
  }

} // namespace Impl

/**
 * \brief VanishingStress material model that enforces stress components to be zero.
 * \ingroup materials
 * \tparam stressIndexPair An array of StressIndexPair representing fixed stress components.
 * \tparam MI The underlying material model.
 */
template <auto stressIndexPair, typename MI>
struct VanishingStress : public Material<VanishingStress<stressIndexPair, MI>>
{
  /**
   * \brief Constructor for VanishingStress.
   * \param mat The underlying material model.
   * \param tol Tolerance for stress reduction.
   */
  explicit VanishingStress(MI mat, typename MI::ScalarType tol = 1e-12)
      : matImpl_{mat},
        tol_{tol} {}

  using Underlying = MI; ///< The underlying material type.

  static constexpr auto fixedPairs        = stressIndexPair;                     ///< Array of fixed stress components.
  static constexpr auto freeVoigtIndices  = createfreeVoigtIndices(fixedPairs);  ///< Free Voigt indices.
  static constexpr auto fixedVoigtIndices = createFixedVoigtIndices(fixedPairs); ///< Fixed Voigt indices.
  static constexpr auto fixedDiagonalVoigtIndicesSize =
      countDiagonalIndices(fixedPairs);                                ///< Number of fixed diagonal indices.
  static constexpr auto freeStrains = freeVoigtIndices.size();         ///< Number of free strains.
  using ScalarType                  = typename Underlying::ScalarType; ///< Scalar type.

  [[nodiscard]] constexpr std::string nameImpl() const noexcept {
    auto matName = matImpl_.name() + "_Vanishing(";
    for (auto p : fixedPairs)
      matName += "(" + std::to_string(p.row) + std::to_string(p.col) + ")";
    matName += ")";
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
  static constexpr double derivativeFactor = 1;                              ///< Derivative factor.

  /**
   * \brief Computes the stored energy for the VanishingStress material.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto [nonOp, Esol] = reduceStress(E);
    return matImpl_.storedEnergyImpl(Esol);
  }

  /**
   * \brief Computes the stresses for the VanishingStress material.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto [nonOp, Esol] = reduceStress(E);
    auto stressesRed         = matImpl_.template stresses<Underlying::strainTag, true>(Esol);

    if constexpr (voigt) {
      return removeCol(stressesRed, fixedVoigtIndices);
    } else
      return fromVoigt(stressesRed, false);
  }

  /**
   * \brief Computes the tangent moduli for the VanishingStress material.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The Green-Lagrangian strain.
   * \return TangentModuli The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& E) const {
    const auto [nonOp, Esol] = reduceStress(E);
    auto C                   = matImpl_.template tangentModuli<Underlying::strainTag, true>(Esol);
    if constexpr (voigt)
      return staticCondensation(C, fixedVoigtIndices);
    else
      return fromVoigt(C);
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam ScalarTypeOther The target scalar type.
   * \return VanishingStress The rebound VanishingStress material.
   */
  template <typename ScalarTypeOther>
  auto rebind() const {
    auto reboundMatImpl = matImpl_.template rebind<ScalarTypeOther>();
    return VanishingStress<stressIndexPair, decltype(reboundMatImpl)>(reboundMatImpl, tol_);
  }

private:
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

  /**
   * \brief Initializes unknown strains based on fixed indices.
   * \tparam Derived The derived type of the input matrix.
   * \param E The input strain matrix.
   */
  template <typename Derived>
  void initUnknownStrains(Eigen::MatrixBase<Derived>& E) const {
    for (size_t i = 0; i < fixedPairs.size(); ++i) {
      ScalarType initialVal = E(fixedPairs[i].row, fixedPairs[i].col);
      if constexpr (strainTag == StrainTags::deformationGradient or strainTag == StrainTags::rightCauchyGreenTensor) {
        if (Dune::FloatCmp::eq(initialVal, ScalarType(0.0)) and (fixedPairs[i].row == fixedPairs[i].col))
          initialVal = ScalarType(1.0);
      }
      if (fixedPairs[i].row != fixedPairs[i].col)
        initialVal = ScalarType(0.0);
      E(fixedPairs[i].row, fixedPairs[i].col) = E(fixedPairs[i].col, fixedPairs[i].row) = initialVal;
    }
  }

  /**
   * \brief Reduces stress components to satisfy the vanishing stress condition.
   * \tparam Derived The derived type of the input matrix.
   * \param Eraw The input strain matrix.
   * \return std::pair<NonLinearOperator, decltype(auto)> The stress reduction result.
   */
  template <typename Derived>
  auto reduceStress(const Eigen::MatrixBase<Derived>& Eraw) const {
    auto E = maybeFromVoigt(Eraw);
    initUnknownStrains(E);

    std::array<size_t, fixedDiagonalVoigtIndicesSize> fixedDiagonalVoigtIndices;
    for (size_t ri = 0; auto i : fixedVoigtIndices) {
      auto indexPair = fromVoigt(i);
      if (indexPair[0] == indexPair[1])
        fixedDiagonalVoigtIndices[ri++] = i;
    }

    auto f = [&](auto&) {
      auto S = matImpl_.template stresses<Underlying::strainTag, true>(E);
      return S(fixedDiagonalVoigtIndices).eval();
    };
    auto df = [&](auto&) {
      auto moduli = (matImpl_.template tangentModuli<Underlying::strainTag, true>(E)).eval();
      return (moduli(fixedDiagonalVoigtIndices, fixedDiagonalVoigtIndices) / Underlying::derivativeFactor).eval();
    };

    auto Er    = E(fixedDiagonalVoigtIndices, fixedDiagonalVoigtIndices).eval().template cast<ScalarType>();
    auto nonOp = Ikarus::NonLinearOperator(functions(f, df), parameter(Er));
    auto nr    = Ikarus::makeNewtonRaphson(
        nonOp, [&](auto& r, auto& A) { return (A.inverse() * r).eval(); },
        [&](auto& /* Ex33 */, auto& ecomps) {
          for (int ri = 0; auto i : fixedDiagonalVoigtIndices) {
            auto indexPair = fromVoigt(i);
            E(indexPair[0], indexPair[1]) += ecomps(ri++);
          }
        });
    nr->setup({.tol = tol_, .maxIter = 100});
    if (!static_cast<bool>(nr->solve()))
      DUNE_THROW(Dune::MathError, "The stress reduction of material " << nameImpl() << " was unsuccessful\n"
                                                                      << "The strains are\n"
                                                                      << E << "\n The stresses are\n"
                                                                      << f(Er));
    return std::make_pair(nonOp, E);
  }

  Underlying matImpl_; ///< The underlying material model.
  double tol_{};       ///< Tolerance for stress reduction.
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
auto makeVanishingStress(MaterialImpl mat, typename MaterialImpl::ScalarType p_tol = 1e-12) {
  return VanishingStress<std::to_array({stressIndexPair...}), MaterialImpl>(mat, p_tol);
}

/**
 * \brief Factory function to create a VanishingStress material for plane stress conditions.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material for plane stress.
 */
template <typename MaterialImpl>
auto planeStress(const MaterialImpl& mat, typename MaterialImpl::ScalarType tol = 1e-8) {
  return makeVanishingStress<Impl::StressIndexPair{2, 1}, Impl::StressIndexPair{2, 0}, Impl::StressIndexPair{2, 2}>(
      mat, tol);
}

/**
 * \brief Factory function to create a VanishingStress material for a shell material with zero normal stress
 * condition.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material for plane stress.
 */
template <typename MaterialImpl>
auto shellMaterial(const MaterialImpl& mat, typename MaterialImpl::ScalarType tol = 1e-8) {
  return makeVanishingStress<Impl::StressIndexPair{2, 2}>(mat, tol);
}

/**
 * \brief Factory function to create a VanishingStress material for a beam material with two zero normal stress
 * condition.
 * \tparam MaterialImpl The underlying material model.
 * \param mat The underlying material model.
 * \param p_tol Tolerance for stress reduction.
 * \return VanishingStress The created VanishingStress material for plane stress.
 */
template <typename MaterialImpl>
auto beamMaterial(const MaterialImpl& mat, typename MaterialImpl::ScalarType tol = 1e-8) {
  return makeVanishingStress<Impl::StressIndexPair{1, 1}, Impl::StressIndexPair{2, 2}>(mat, tol);
}
} // namespace Ikarus
