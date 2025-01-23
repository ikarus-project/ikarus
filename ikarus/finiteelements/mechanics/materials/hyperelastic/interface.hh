// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Implementation of the Hyperelastic material model.
 * \ingroup  materials
 */

#pragma once

#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/volumetric/volumetricfunctions.hh>
#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of a general Hyperelastic Material material model.
 * \details \f$\Psi(\BC) = \hat{\Psi}(\la_1, \la_2, \la_3) + U(J)\f$ with \f$\hat{\Psi}\f$ being the
 *deviatoric part of the strain energy function and \f$ U(J) \f$ being the volumetric part. After calling the underlying
 *deviatoric and volumetric function, the transformation to cartesian coordinate system is implemented in this
 *interface.
 *
 * \ingroup materials
 */
template <typename DEV, typename VOL = NoVolumetricPart>
requires(std::same_as<typename DEV::ScalarType, typename VOL::ScalarType>)
struct Hyperelastic : public Material<Hyperelastic<DEV, VOL>>
{
  static_assert(std::is_same_v<DEV, Deviatoric<typename DEV::DeviatoricFunction>>);
  static_assert(std::is_same_v<VOL, Volumetric<typename VOL::VolumetricFunction>>);

  using ScalarType                        = typename DEV::ScalarType;
  static constexpr bool hasVolumetricPart = not std::same_as<VOL, NoVolumetricPart>;
  static constexpr bool isAutoDiff        = Concepts::AutodiffScalar<ScalarType>;

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
  static constexpr bool energyAcceptsVoigt     = false;
  static constexpr bool stressToVoigt          = false;
  static constexpr bool stressAcceptsVoigt     = false;
  static constexpr bool moduliToVoigt          = false;
  static constexpr bool moduliAcceptsVoigt     = false;
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
        vol_(NoVolumetricPart{MaterialParametersVOL{}, {}}) {}

  Hyperelastic(const DEV& dev, const VOL& vol)
      : dev_(dev),
        vol_(vol) {}

  /** \brief Returns the deviatoric function. */
  const DEV& deviatoricFunction() const { return dev_; }

  /** \brief Returns the volumetric function. */
  const VOL& volumetricFunction() const { return vol_; }

  /**
   * \brief Returns the material parameters stored in the deviatoric part of the material.
   */
  const MaterialParameters materialParametersImpl() const {
    if constexpr (hasVolumetricPart)
      return {dev_.materialParameters(), vol_.materialParameter()};
    else
      return dev_.materialParameters();
  }

  /**
   * \brief Computes the total stored energy in the Hyperelastic material model.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return ScalarType The stored energy.
   */
  template <typename Derived>
  ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);

    const auto lambdas = principalStretches(C, Eigen::EigenvaluesOnly).first;

    // Workaround to avoid the usage of degenerated principal stretches while using AutoDiff
    ScalarType J = isAutoDiff ? sqrt(C.derived().eval().determinant()) : detF(lambdas);

    const auto devEnergy =
        isAutoDiff ? deviatoricEnergy(C)
                   : dev_.storedEnergy(lambdas); // avoid duplicate computation of principal stretches for performance

    return devEnergy + vol_.storedEnergy(J);
  }

  /**
   * \brief Computes the stresses in the Hyperelastic material model.
   * \tparam voigt A boolean indicating whether to return stresses in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return StressMatrix The stresses.
   */
  template <bool voigt, typename Derived>
  StressMatrix stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      const auto [lambdas, N] = principalStretches(C);

      // Workaround to avoid the usage of degenerated principal stretches while using AutoDiff
      ScalarType J = isAutoDiff ? sqrt(C.derived().eval().determinant()) : detF(lambdas);

      const auto Sdev = isAutoDiff ? deviatoricStress(C)
                                   : transformDeviatoricStresses(
                                         dev_.stresses(lambdas),
                                         N); // avoid duplicate computation of principal stretches for performance
      const auto Svol = transformVolumetricStresses(vol_.firstDerivative(J), C, J);

      return Sdev + Svol;
    } else
      static_assert(voigt == false, "Hyperelastic does not support returning stresses in Voigt notation");
  }

  /**
   * \brief Computes the tangent moduli in the Hyperelastic material model.
   * \tparam voigt A boolean indicating whether to return tangent moduli in Voigt notation.
   * \tparam Derived The derived type of the input matrix.
   * \param C The right Cauchy-Green tensor.
   * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> The tangent moduli.
   */
  template <bool voigt, typename Derived>
  MaterialTensor tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      const auto [lambdas, N] = principalStretches(C);

      // Workaround to avoid the usage of degenerated principal stretches while using AutoDiff
      ScalarType J = isAutoDiff ? sqrt(C.derived().eval().determinant()) : detF(lambdas);

      const auto moduliDev = transformDeviatoricTangentModuli(dev_.tangentModuli(lambdas), N);
      const auto moduliVol = transformVolumetricTangentModuli(vol_.firstDerivative(J), vol_.secondDerivative(J), C, J);

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

  inline static constexpr auto dimensionRange() { return Dune::range(dim); }

  template <typename ST>
  Eigen::Matrix<ST, 3, 3> transformDeviatoricStresses(const Eigen::Matrix<ST, 3, 1>& principalStress,
                                                      const Eigen::Matrix<ST, 3, 3>& N) const {
    return (N * principalStress.asDiagonal() * N.transpose()).eval();
  }

  template <typename ST>
  Eigen::Matrix<ST, 3, 3> transformVolumetricStresses(ST Uprime, const auto& C, ST J) const {
    return J * Uprime * C.inverse();
  }

  template <typename ST>
  auto transformDeviatoricTangentModuli(const Eigen::TensorFixedSize<ST, Eigen::Sizes<3, 3, 3, 3>>& L,
                                        const Eigen::Matrix<ST, 3, 3>& N) const {
    Eigen::TensorFixedSize<ST, Eigen::Sizes<3, 3, 3, 3>> moduli{};
    moduli.setZero();

    for (const auto i : dimensionRange())
      for (const auto k : dimensionRange()) {
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

  template <typename ST>
  auto transformVolumetricTangentModuli(const ST& Uprime, const ST& Uprimeprime, const auto& C, ST J) const {
    const auto invC  = C.inverse().eval();
    const auto CTinv = tensorView(invC, std::array<Eigen::Index, 2>({3, 3}));

    Eigen::TensorFixedSize<ST, Eigen::Sizes<3, 3, 3, 3>> moduli =
        (J * ((Uprime + J * Uprimeprime) * dyadic(CTinv, CTinv) -
              (2 * Uprime * symTwoSlots(fourthOrderIKJL(invC, invC), {2, 3}))))
            .eval();

    return moduli;
  }

  template <typename Derived>
  auto principalStretches(const Eigen::MatrixBase<Derived>& Craw, int options = Eigen::ComputeEigenvectors) const {
    auto C = Impl::maybeFromVoigt(Craw);
    return Impl::principalStretches<typename Derived::Scalar>(C, options);
  }

  template <typename ST>
  auto detF(const Eigen::Array<ST, 3, 1>& lambda) const -> ST {
    if constexpr (hasVolumetricPart) {
      const auto detC = Impl::determinantFromPrincipalValues<ST>(lambda);
      Impl::checkPositiveDet(detC);
      return detC;
    }
    return 0.0;
  }

  // Unpack the derivatives from the result of an @ref eval call into an array.
  template <typename D, typename F>
  auto forEach(const Eigen::MatrixBase<D>& result, F&& f) const {
    auto& r = result.derived();
    return r.unaryExpr(f);
  }

  /** \brief A helper function to compute the deviatoric part of the energy.
   *
   * \details While using AutoDiff, if the eigenvalues (principal stretches) of C are degenerated, the derivative can
   * have certain singularities. This is circumvented here by explicitly updating the derivatives.
   *
   */
  template <typename Derived>
  auto deviatoricEnergy(const Eigen::MatrixBase<Derived>& C) const {
    if constexpr (not Concepts::AutodiffScalar<typename Derived::Scalar>) {
      auto [lambdas, N] = principalStretches(C);
      return dev_.storedEnergy(lambdas);
    } else if constexpr (std::is_same_v<ScalarType, autodiff::dual>) {
      autodiff::dual e;
      auto Cvec     = toVoigt(C.derived());
      auto realCVec = autodiff::derivative<0>(Cvec);
      auto dualCVec = autodiff::derivative<1>(Cvec);

      auto [lambdas, N] = principalStretches(realCVec);

      auto realDev = dev_.template rebind<double>();

      e.val = realDev.storedEnergy(lambdas);
      e.grad =
          (transformDeviatoricStresses(realDev.stresses(lambdas), N).transpose() / 2 * fromVoigt(dualCVec)).trace();
      return e;
    } else if constexpr (std::is_same_v<ScalarType, autodiff::dual2nd>) {
      autodiff::dual2nd e;
      auto Cvec           = toVoigt(C.derived());
      const auto realCVec = derivative<0>(Cvec);
      const auto dualC    = fromVoigt(forEach(Cvec, [](auto& v) { return v.grad.val; }).eval());
      const auto dualC2   = fromVoigt(forEach(Cvec, [](auto& v) { return v.val.grad; }).eval());
      auto [lambdas, N]   = principalStretches(realCVec);

      auto realDev = dev_.template rebind<double>();
      e.val        = realDev.storedEnergy(lambdas);

      e.grad.val         = (transformDeviatoricStresses(realDev.stresses(lambdas), N).transpose() / 2 * dualC).trace();
      e.val.grad         = e.grad.val;
      const auto Cmoduli = transformDeviatoricTangentModuli(realDev.tangentModuli(lambdas), N);

      Eigen::array<Eigen::IndexPair<Eigen::Index>, 2> double_contraction  = {Eigen::IndexPair<Eigen::Index>(2, 0),
                                                                             Eigen::IndexPair<Eigen::Index>(3, 1)};
      Eigen::array<Eigen::IndexPair<Eigen::Index>, 2> double_contraction2 = {Eigen::IndexPair<Eigen::Index>(0, 0),
                                                                             Eigen::IndexPair<Eigen::Index>(1, 1)};
      const auto tCdual                  = tensorView(dualC, std::array<Eigen::Index, 2>({3, 3}));
      const auto tCdualT                 = tensorView(dualC2, std::array<Eigen::Index, 2>({3, 3}));
      const auto prod                    = Cmoduli.contract(tCdual, double_contraction);
      const Eigen::Tensor<double, 0> res = tCdualT.contract(prod, double_contraction2);
      e.grad.grad                        = res(0) / 4.0; // extracting value of zero order tensor

      return e;
    } else
      static_assert(Dune::AlwaysFalse<Derived>::value, "No fitting ScalarType.");
  }

  /** \brief A helper function to compute the deviatoric part of the stresses.
   *
   * \details While using AutoDiff, if the eigenvalues (principal stretches) of C are degenerated, the derivative can
   * have certain singularities. This is circumvented here by explicitly updating the derivatives.
   *
   */
  template <typename Derived>
  auto deviatoricStress(const Eigen::MatrixBase<Derived>& C) const {
    if constexpr (not Concepts::AutodiffScalar<typename Derived::Scalar>) {
      auto [lambdas, N] = principalStretches(C);
      return transformDeviatoricStresses(dev_.stresses(lambdas), N);
    } else if constexpr (std::is_same_v<ScalarType, autodiff::dual>) {
      constexpr int nVoigtIndices = 6;
      Eigen::Vector<autodiff::dual, nVoigtIndices> g;
      auto Cvec           = toVoigt(C.derived());
      const auto realCVec = derivative<0>(Cvec);
      auto realC          = fromVoigt(realCVec);
      auto dualC          = fromVoigt(forEach(Cvec, [](const auto& v) { return v.grad; }).eval());
      auto [lambdas, N]   = principalStretches(realC);

      auto realDev       = dev_.template rebind<double>();
      const auto Cmoduli = toVoigt(transformDeviatoricTangentModuli(realDev.tangentModuli(lambdas), N));
      for (int i = 0; i < nVoigtIndices; ++i) {
        Eigen::Vector<double, nVoigtIndices> contraction = Cmoduli * toVoigt(dualC);
        contraction.topRows<3>() /= 2.0;
        g[i].val  = toVoigt(transformDeviatoricStresses(realDev.stresses(lambdas), N))[i];
        g[i].grad = contraction[i];
      }

      return fromVoigt(g);
    } else
      static_assert(Dune::AlwaysFalse<Derived>::value, "No fitting ScalarType.");
  }
};

} // namespace Ikarus::Materials
