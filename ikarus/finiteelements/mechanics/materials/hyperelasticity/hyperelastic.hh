// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file hyperelastic.hh
 * \brief Implementation of the Hyperelastic material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/concepts.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelasticity/volumetric.hh>
#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of a general Hyperelastic Material material model.
 * \details $\Psi(\BC) = \hat{\Psi}(\lambda_1 + \lambda_2 + \lambda_3) + U(J)$ with $\hat{\Psi}$ being the deviatoric
 * part of the strain energy function and U(J) the volumetric part.
 * \ingroup materials
 */
template <typename DEV, typename VOL = NoVolumetricPart>
requires(std::same_as<typename DEV::ScalarType, typename VOL::ScalarType>)
struct Hyperelastic : public Material<Hyperelastic<DEV, VOL>>
{
  // Checking concepts here results in better compiler error messages because of CRTP
  static_assert(Concepts::DeviatoricPart<DEV>);
  static_assert(Concepts::VolumetricPart<VOL>);

  using ScalarType                        = typename DEV::ScalarType;
  static constexpr bool hasVolumetricPart = not std::same_as<VOL, NoVolumetricPart>;

  static constexpr int dim = 3;
  using StrainMatrix       = Eigen::Matrix<ScalarType, dim, dim>;
  using StressMatrix       = StrainMatrix;
  using MaterialTensor     = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>>;

  using MaterialParametersDEV = typename DEV::MaterialParameters;
  using MaterialParametersVOL = typename VOL::MaterialParameter;
  using MaterialParameters =
      std::conditional_t<hasVolumetricPart, std::pair<MaterialParametersDEV, MaterialParametersVOL>,
                         MaterialParametersDEV>;

  static constexpr auto strainTag              = StrainTags::rightCauchyGreenTensor;
  static constexpr auto stressTag              = StressTags::PK2;
  static constexpr auto tangentModuliTag       = TangentModuliTags::Material;
  static constexpr bool energyAcceptsVoigt     = true;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = true;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = true;
  static constexpr double derivativeFactorImpl = 2;

  [[nodiscard]] constexpr static std::string nameImpl() noexcept {
    if constexpr (hasVolumetricPart)
      return "Hyperelastic (" + DEV::name() + ", " + VOL::name() + ")";
    else
      return "Hyperelastic (" + DEV::name() + ")";
  }

  explicit Hyperelastic(const DEV& dev)
  requires(not hasVolumetricPart)
      : dev_{dev},
        vol_(NoVolumetricPart{BulkModulus{}, {}}) {}

  Hyperelastic(const DEV& dev, const VOL& vol)
      : dev_(dev),
        vol_(vol) {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return dev_.materialParameter_; }

  /**
   * \brief Computes the stored energy in the Neo-Hookean material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);

    auto lambdas = principalStretches(C, Eigen::EigenvaluesOnly).first;
    auto J       = detF(C);

    return dev_.storedEnergy(lambdas) + vol_.storedEnergy(J);
  }

  /**
   * \brief Computes the stresses in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  StressMatrix stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      auto [lambdas, N] = principalStretches(C);
      auto J            = detF(C);

      auto Sdev = transformDeviatoricStresses(dev_.stresses(lambdas), N);
      auto Svol = transformVolumetricStresses(vol_.firstDerivative(J), C, J);

      return Sdev + Svol;
    } else
      static_assert(voigt == false, "Hyperelastic does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Neo-Hookean material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  MaterialTensor tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      auto [lambdas, N] = principalStretches(C);
      auto J            = detF(C);

      auto moduliDev = transformDeviatoricTangentModuli(dev_.tangentModuli(lambdas), N);
      auto moduliVol = transformVolumetricTangentModuli(vol_.firstDerivative(J), vol_.secondDerivative(J), C, J);

      return moduliDev + moduliVol;
    } else
      static_assert(voigt == false, "Hyperelastic does not support returning tangent moduli in Voigt notation");
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return Hyperelastic<ScalarTypeOther> The rebound Hyperelastic material.
   */
  template <typename STO>
  auto rebind() const {
    auto reboundDEV = dev_.template rebind<STO>();
    auto reboundVOL = vol_.template rebind<STO>();
    return Hyperelastic<decltype(reboundDEV), decltype(reboundVOL)>(reboundDEV, reboundVOL);
  }

private:
  DEV dev_;
  VOL vol_;

  StressMatrix transformDeviatoricStresses(const typename DEV::StressMatrix& principalStress, const auto& N) const {
    auto S = StressMatrix::Zero().eval();
    for (auto i : Dune::range(3))
      S += principalStress[i] * dyadic<ScalarType, 3, false>(N.col(i).eval(), N.col(i).eval());

    return S;
  }

  StressMatrix transformVolumetricStresses(const ScalarType& Uprime, const auto& C, ScalarType J) const {
    return J * Uprime * C.inverse();
  }

  MaterialTensor transformDeviatoricTangentModuli(const typename DEV::MaterialTensor& L, const auto& N) const {
    MaterialTensor moduli{};
    moduli.setZero();

    for (int i = 0; i < dim; ++i)
      for (int k = 0; k < dim; ++k) {
        // First term: L[i, i, k, k] * ((N[i] ⊗ N[i]) ⊗ (N[k] ⊗ N[k]))
        auto NiNi = dyadic(N.col(i).eval(), N.col(i).eval());
        auto NkNk = dyadic(N.col(k).eval(), N.col(k).eval());

        moduli += L(i, i, k, k) * dyadic(NiNi, NkNk);

        // Second term (only if i != k): L[i, k, i, k] * (N[i] ⊗ N[k] ⊗ (N[i] ⊗ N[k] + N[k] ⊗ N[i]))
        if (i != k) {
          auto NiNk = dyadic(N.col(i).eval(), N.col(k).eval());
          auto NkNi = dyadic(N.col(k).eval(), N.col(i).eval());

          moduli += L(i, k, i, k) * dyadic(NiNk, NiNk + NkNi);
        }
      }

    return moduli;
  }

  MaterialTensor transformVolumetricTangentModuli(const ScalarType& Uprime, const ScalarType& Uprimeprime,
                                                  const auto& C, ScalarType J) const {
    const auto invC    = C.inverse().eval();
    const auto CTinv   = tensorView(invC, std::array<Eigen::Index, 2>({3, 3}));
    const auto CinvDya = dyadic(CTinv, CTinv);
    const auto CinvT23 = symTwoSlots(fourthOrderIKJL(invC, invC), {2, 3});

    Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> moduli =
        (J * ((Uprime + J * Uprimeprime) * CinvDya - 2 * Uprime * CinvT23)).eval();

    return moduli;
  }

  /**
   * \brief Calculates the principal stretches of the input strain matrix C
   *
   * \tparam Derived The derived type of the input matrix
   * \param Craw the input strain matrix
   * \param options should be either `Eigen::ComputeEigenvectors` or `Eigen::EigenvaluesOnly`
   * \return auto pair of principalstretches and corresponding eigenvectors (if `Eigen::EigenvaluesOnly` the
   * eigenvectors is a zero matrix )
   */
  template <typename Derived>
  auto principalStretches(const Eigen::MatrixBase<Derived>& Craw, int options = Eigen::ComputeEigenvectors) const {
    StrainMatrix C = Impl::maybeFromVoigt(Craw);
    Eigen::SelfAdjointEigenSolver<StrainMatrix> eigensolver{};

    // For AD we use the direct method which uses a closed form calculation, but has less accuracy
    if constexpr (Concepts::AutodiffScalar<ScalarType>)
      eigensolver.computeDirect(C, options);
    else
      eigensolver.compute(C, options);

    auto& eigenvalues  = eigensolver.eigenvalues();
    auto& eigenvectors = options == Eigen::ComputeEigenvectors ? eigensolver.eigenvectors() : StrainMatrix::Zero();

    auto principalStretches = eigenvalues.array().sqrt().eval();
    return std::make_pair(principalStretches, eigenvectors);
  }

  template <typename Derived>
  auto detF(const Eigen::MatrixBase<Derived>& Craw) const -> typename VOL::JType {
    StrainMatrix C = Impl::maybeFromVoigt(Craw);

    if constexpr (hasVolumetricPart) {
      const auto detC = C.determinant();
      Impl::checkPositiveDet(detC);

      return sqrt(detC);
    }
    return 0.0;
  }
};

} // namespace Ikarus::Materials