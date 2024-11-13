// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file Muesli.hh
 * \brief Implementation of the Muesli material model.
 * \ingroup  materials
 */

#pragma once

#include <muesli/muesli.h>

#include <Eigen/Eigen>

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslihelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

template <typename FM>
requires(std::is_base_of_v<muesli::finiteStrainMaterial, FM>)
struct MuesliFinite : public Material<MuesliFinite<FM>>
{
  using MaterialModel                 = FM;
  using ScalarType                    = double;
  static constexpr int worldDimension = 3;
  using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
  using StressMatrix                  = StrainMatrix;
  using MaterialParameters            = Muesli::MaterialProperties;

  static constexpr auto strainTag              = StrainTags::deformationGradient;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
  static constexpr double derivativeFactorImpl = 1;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept { return "Muesli: " + Dune::className<FM>(); }

  /**
   * \brief Constructor for MuesliT.
   * \param mpt The Lame's parameters (first parameter and shear modulus).
   */
  explicit MuesliFinite(const MaterialParameters& mpt)
      : materialParameter_{mpt},
        material_{Dune::className<FM>(), mpt},
        mp_{material_.createMaterialPoint()} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameter_; }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!Concepts::EigenVector<Derived>) {
      istensor strain = istensor(C(0, 0), C(1, 1), C(2, 2), C(1, 2), C(2, 0), C(0, 1));

      mp_->updateCurrentState(0.0, strain);
      return mp_->storedEnergy();

    } else
      static_assert(!Concepts::EigenVector<Derived>,
                    "Muesli energy can only be called with a matrix and not a vector in Voigt notation");
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
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
        istensor strain = istensor(C(0, 0), C(1, 1), C(2, 2), C(1, 2), C(2, 0), C(0, 1));

        mp_->updateCurrentState(0.0, strain);

        auto stress = istensor();
        mp_->secondPiolaKirchhoffStress(stress);

        auto S = Eigen::Matrix<double, 3, 3>{};
        for (auto i : Dune::range(3))
          for (auto j : Dune::range(3))
            S(i, j) = stress(i, j);
        return S;

      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "Muesli can only be called with a matrix and not a vector in Voigt notation");
    } else
      static_assert(voigt == false, "Muesli does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      istensor strain = istensor(C(0, 0), C(1, 1), C(2, 2), C(1, 2), C(2, 0), C(0, 1));

      mp_->updateCurrentState(0.0, strain);

      Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> moduli{};
      moduli.setZero();

      auto tangent = itensor4();
      mp_->convectedTangent(tangent);

      for (auto i : Dune::range(3))
        for (auto j : Dune::range(3))
          for (auto k : Dune::range(3))
            for (auto l : Dune::range(3))
              moduli(i, j, k, l) = tangent(i, j, k, l);

      return moduli;
    } else
      static_assert(voigt == false, "Muesli does not support returning tangent moduli in Voigt notation");
  }

  auto& material() const { return material_; }

  ~MuesliFinite() { delete mp_; }

  MuesliFinite(const MuesliFinite& other)
      : materialParameter_{other.materialParameter_},
        material_{Dune::className<FM>(), materialParameter_},
        mp_{material_.createMaterialPoint()} {}

private:
  MaterialParameters materialParameter_;
  MaterialModel material_;
  muesli::finiteStrainMP* mp_;
};

} // namespace Ikarus::Materials
