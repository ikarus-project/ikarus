// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file Muesli.hh
 * \brief Implementation of the Muesli material model.
 * \ingroup materials
 */

#pragma once

#if ENABLE_MUESLI

  #include <muesli/muesli.h>

  #include <Eigen/Eigen>

  #include <ikarus/finiteelements/mechanics/materials/interface.hh>
  #include <ikarus/finiteelements/mechanics/materials/muesli/mueslihelpers.hh>
  #include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Wrapper class for small strain materials from the muesli library. It can be templated with all materials
 * derived from muesli::elasticIsotropicMaterial. It models the Ikarus material interface.
 *
 * \tparam SM muesli material model implementation
 */
template <typename SM = muesli::elasticIsotropicMaterial>
requires(std::is_base_of_v<muesli::smallStrainMaterial, SM>)
struct SmallStrain : public Material<SmallStrain<SM>>
{
  using MaterialModel      = SM;
  using ScalarType         = double;
  static constexpr int dim = 3;
  using StrainMatrix       = Eigen::Matrix<ScalarType, dim, dim>;
  using StressMatrix       = StrainMatrix;
  using MaterialTensor     = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>>;
  using MaterialParameters = muesli::materialProperties;

  static constexpr auto strainTag              = StrainTags::linear;
  static constexpr auto stressTag              = StressTags::linear;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 1;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "SmallStrain: " + materialName<SM>(); }

  /**
   * \brief Constructor for small strain muesli materials (only activated for isotropic linear elasticity).
   * \param mpt Arbitrary Material Parameter tuple defined in \file physicshelper.hh.
   */
  template <Concepts::MPTuple MPT>
  requires(std::same_as<MaterialModel, muesli::elasticIsotropicMaterial>)
  explicit SmallStrain(const MPT& mpt)
      : materialParameter_{propertiesFromIkarusMaterialParameters(mpt)},
        material_{Dune::className<SM>(), materialParameter_},
        mp_{material_.createMaterialPoint()} {}

  /**
   * \brief Constructor for small strain muesli materials.
   * \param mpt Muesli materialproperties.
   */
  explicit SmallStrain(const MaterialParameters& mpt)
      : materialParameter_{mpt},
        material_{Dune::className<SM>(), mpt},
        mp_{material_.createMaterialPoint()} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  /**
   * \brief Computes the stored energy in the Muesli small strain material model.
   * \tparam Derived The derived type of the input matrix.
   * \param E The linear strain tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!Concepts::EigenVector<Derived>) {
      updateState(E);
      return mp_->storedEnergy();
    } else
      static_assert(!Concepts::EigenVector<Derived>,
                    "MuesliSmallStrain energy can only be called with a matrix and not a vector in Voigt notation");
  }

  /**
   * \brief Computes the stresses in the Muesli small strain material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The linear strain tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  auto stressesImpl(const Eigen::MatrixBase<Derived>& E) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        updateState(E);
        istensor stress;
        mp_->stress(stress);
        return toEigenMatrix(stress);
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "MuesliSmallStrain can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "MuesliSmallStrain does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Muesli small strain material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param E The linear strain tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& E) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      if constexpr (!Concepts::EigenVector<Derived>) {
        updateState(E);
        itensor4 tangentModuli;
        mp_->tangentTensor(tangentModuli);
        return toEigenTensor(tangentModuli);
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "MuesliSmallStrain can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "MuesliSmallStrain does not support returning tangent moduli in Voigt notation");
  }

  /**
   * \brief Returns the underlying muesli material implementation
   * \return auto& reference to the muesli material
   */
  auto& material() const { return material_; }

  /**
   * \brief Returns the underlying muesli material point implementation.
   * \return auto& reference to the muesli material point.
   */
  auto& materialPoint() const { return mp_; }

  SmallStrain(const SmallStrain& other)
      : materialParameter_{other.materialParameter_},
        material_{Dune::className<SM>(), materialParameter_},
        mp_{material_.createMaterialPoint()} {}

private:
  MaterialParameters materialParameter_;
  MaterialModel material_;
  std::unique_ptr<muesli::smallStrainMP> mp_;

  template <typename Derived>
  void updateState(const Eigen::MatrixBase<Derived>& E) const {
    mp_->updateCurrentState(0.0, toistensor(E.derived()));
  }
};

} // namespace Ikarus::Materials

#else
  #error Muesli materials depends on the Muesli library, which is not included
#endif