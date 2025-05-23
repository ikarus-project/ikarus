// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file Muesli.hh
 * \brief Implementation of the Muesli material model for finite strains.
 * \ingroup materials
 */

#pragma once

#if ENABLE_MUESLI

  #include <muesli/muesli.h>

  #include <Eigen/Eigen>

  #include <ikarus/finiteelements/mechanics/materials/interface.hh>
  #include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
  #include <ikarus/finiteelements/mechanics/materials/muesli/mueslihelpers.hh>
  #include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Wrapper class for finite strain materials from the muesli library. It can be templated with all materials
 * derived from muesli::finiteStrainMaterial. It adheres to the Ikarus material interface. See
 * \file ikarus/finiteelements/mechanics/materials/interface.hh.
 *
 * \tparam FM muesli material model implementation
 */
template <typename FM>
requires(std::is_base_of_v<muesli::finiteStrainMaterial, FM>)
struct FiniteStrain : public Material<FiniteStrain<FM>>
{
  using MaterialModel      = FM;
  using ScalarType         = double;
  static constexpr int dim = 3;
  using StrainMatrix       = Eigen::Matrix<ScalarType, dim, dim>;
  using StressMatrix       = StrainMatrix;
  using MaterialTensor     = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>>;
  using MaterialParameters = muesli::materialProperties;

  static constexpr auto strainTag              = StrainTags::rightCauchyGreenTensor;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 2;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "FiniteStrain: " + materialName<FM>(); }

  /**
   * \brief Constructor for finite strain muesli materials.
   * \param mpt Muesli materialproperties.
   */
  explicit FiniteStrain(const MaterialParameters& mpt)
      : materialParameter_{mpt},
        material_{Dune::className<FM>(), mpt},
        mp_{material_.createMaterialPoint()} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  /**
   * \brief Computes the stored energy in the Muesli finite strain material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!Concepts::EigenVector<Derived>) {
      updateState(C);
      return mp_->storedEnergy();
    } else
      static_assert(!Concepts::EigenVector<Derived>,
                    "MuesliFiniteStrain energy can only be called with a matrix and not a vector in Voigt notation");
  }

  /**
   * \brief Computes the stresses in the Muesli finite strain material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        updateState(C);
        istensor stress;
        mp_->secondPiolaKirchhoffStress(stress);
        return toEigenMatrix(stress).eval();
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "MuesliFiniteStrain can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "MuesliFiniteStrain does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Muesli finite strain material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        updateState(C);
        itensor4 tangentModuli;

        mp_->convectedTangent(tangentModuli);
        return toEigenTensor(tangentModuli);
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "MuesliFiniteStrain can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "MuesliFiniteStrain does not support returning tangent moduli in Voigt notation");
  }

  /**
   * \brief Returns the underlying muesli material implementation.
   * \return auto& reference to the muesli material.
   */
  auto& material() const { return material_; }

  /**
   * \brief Returns the underlying muesli material point implementation.
   * \return auto& reference to the muesli material point.
   */
  auto& materialPoint() const { return mp_; }

  FiniteStrain(const FiniteStrain& other)
      : materialParameter_{other.materialParameter_},
        material_{Dune::className<FM>(), materialParameter_},
        mp_{material_.createMaterialPoint()} {}

private:
  MaterialParameters materialParameter_;
  MaterialModel material_;
  std::unique_ptr<muesli::finiteStrainMP> mp_;

  template <typename Derived>
  void updateState(const Eigen::MatrixBase<Derived>& C) const {
    Impl::checkPositiveOrAbort(C.determinant());
    auto F = toitensor(transformStrain<strainTag, StrainTags::deformationGradient>(C));
    mp_->updateCurrentState(0.0, F);
  }
};

} // namespace Ikarus::Materials
#else
  #error Muesli materials depends on the Muesli library, which is not included
#endif