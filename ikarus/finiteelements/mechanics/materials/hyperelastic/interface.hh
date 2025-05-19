// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
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
 * deviatoric part of the strain energy function and \f$ U(J) \f$ being the volumetric part. After calling the
 * underlying deviatoric and volumetric function, the transformation to cartesian coordinate system is implemented in
 * this interface.
 *
 * \ingroup materials
 */
template <typename DEV, typename VOL = NoVolumetricPart>
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
  using MaterialTensor     = Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>>;

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
        vol_(VOL{MaterialParametersVOL{}, typename VOL::VolumetricFunction{}}) {}

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
    auto J             = detF(C, lambdas);

    return deviatoricEnergy(C, lambdas) + vol_.storedEnergy(J);
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
      auto J                  = detF(C, lambdas);

      const auto Sdev = deviatoricStress(C, lambdas, N);
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
   * \return MaterialTensor The tangent moduli.
   */
  template <bool voigt, typename Derived>
  MaterialTensor tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
    static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
    if constexpr (!voigt) {
      const auto [lambdas, N] = principalStretches(C);
      auto J                  = detF(C, lambdas);

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
    return Hyperelastic<decltype(reboundDEV), VOL>(reboundDEV, vol_);
  }

private:
  DEV dev_;
  VOL vol_;

  inline static constexpr auto dimensionRange() { return Dune::range(dim); }

  template <typename ST>
  Eigen::Matrix<ST, dim, dim> transformDeviatoricStresses(const Eigen::Vector<ST, dim>& principalStress,
                                                          const Eigen::Matrix<ST, dim, dim>& N) const {
    return (N * principalStress.asDiagonal() * N.transpose()).eval();
  }

  template <typename ST>
  Eigen::Matrix<ST, dim, dim> transformVolumetricStresses(ST Uprime, const auto& C, ST J) const {
    return J * Uprime * C.inverse();
  }

  template <typename ST>
  auto transformDeviatoricTangentModuli(const Eigen::TensorFixedSize<ST, Eigen::Sizes<dim, dim, dim, dim>>& L,
                                        const Eigen::Matrix<ST, dim, dim>& N) const {
    Eigen::TensorFixedSize<ST, Eigen::Sizes<dim, dim, dim, dim>> moduli{};
    moduli.setZero();
    auto indexList = std::array<Eigen::Index, 2>({dim, dim});

    for (const auto i : dimensionRange())
      for (const auto k : dimensionRange()) {
        // First term: L[i, i, k, k] * ((N[i] ⊗ N[i]) ⊗ (N[k] ⊗ N[k]))
        auto NiNi = tensorView(dyadic(N.col(i).eval(), N.col(i).eval()), indexList);
        auto NkNk = tensorView(dyadic(N.col(k).eval(), N.col(k).eval()), indexList);

        moduli += L(i, i, k, k) * dyadic(NiNi, NkNk);

        // Second term (only if i != k): L[i, k, i, k] * (N[i] ⊗ N[k] ⊗ (N[i] ⊗ N[k] + N[k] ⊗ N[i]))
        if (i != k) {
          auto NiNk = tensorView(dyadic(N.col(i).eval(), N.col(k).eval()), indexList);
          auto NkNi = tensorView(dyadic(N.col(k).eval(), N.col(i).eval()), indexList);

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

  /**
   * \brief A helper function to determine the determinant of the deformation gradient F.
   *
   * \details According to https://eigen.tuxfamily.org/dox/Tridiagonalization_8h_source.html (line 477), which is used
   * while computing the eigenvalues to obtain the prinipal stretches, if the (2,0) entry of a matrix is zero, an early
   * exit is performed. Even though this doesn't affect the results in the case of doubles, using an autodiff type leads
   * to problems as information is lost. Alternatively, using `computeDirect()` to compute the eigenvalues leads to NaNs
   * as sqrt of zero is being taken. Therefore, if autodiff is used, the determinant is computed explicitly from the
   * right Cauchy-Green tensor C.
   */
  template <typename Derived, typename ST>
  ST detF(const Eigen::MatrixBase<Derived>& C, const Eigen::Vector<ST, 3>& lambda) const {
    if constexpr (isAutoDiff) {
      const auto detC = sqrt(C.derived().eval().determinant());
      return detC;
    } else {
      const auto detC = Impl::determinantFromPrincipalValues(lambda);
      Impl::checkPositiveOrAbort(detC);
      return detC;
    }
  }

  /** \brief A helper function to compute the deviatoric part of the energy.
   *
   * \details While using AutoDiff, if the eigenvalues (principal stretches) of C are degenerated, the derivative can
   * have certain singularities. This is circumvented here by explicitly updating the derivatives.
   *
   */
  template <typename Derived, typename ST>
  requires(std::same_as<typename Derived::Scalar, ST>)
  auto deviatoricEnergy(const Eigen::MatrixBase<Derived>& C, const Eigen::Vector<ST, 3>& lambdasST) const {
    if constexpr (not Concepts::AutodiffScalar<ST>) {
      return dev_.storedEnergy(lambdasST);
    } else if constexpr (std::is_same_v<ST, autodiff::dual>) {
      autodiff::dual e;
      const auto Cvec     = toVoigt(C.derived());
      const auto realCVec = autodiff::derivative<0>(Cvec);
      const auto dualCVec = autodiff::derivative<1>(Cvec);

      const auto [lambdas, N] = principalStretches(realCVec);

      e.val  = dev_.storedEnergy(lambdas);
      e.grad = (transformDeviatoricStresses(dev_.stresses(lambdas), N).transpose() / 2 * fromVoigt(dualCVec)).trace();
      return e;
    } else if constexpr (std::is_same_v<ST, autodiff::dual2nd>) {
      autodiff::dual2nd e;
      const auto Cvec         = toVoigt(C.derived());
      const auto realCVec     = derivative<0>(Cvec);
      const auto dualC        = fromVoigt(Cvec.unaryExpr([](auto& v) { return v.grad.val; }).eval());
      const auto dualC2       = fromVoigt(Cvec.unaryExpr([](auto& v) { return v.val.grad; }).eval());
      const auto [lambdas, N] = principalStretches(realCVec);

      e.val      = dev_.storedEnergy(lambdas);
      e.grad.val = (transformDeviatoricStresses(dev_.stresses(lambdas), N).transpose() / 2 * dualC).trace();
      e.val.grad = e.grad.val;

      const auto Cmoduli = transformDeviatoricTangentModuli(dev_.tangentModuli(lambdas), N);
      const Eigen::array<Eigen::IndexPair<Eigen::Index>, 2> double_contraction  = {Eigen::IndexPair<Eigen::Index>(2, 0),
                                                                                   Eigen::IndexPair<Eigen::Index>(3, 1)};
      const Eigen::array<Eigen::IndexPair<Eigen::Index>, 2> double_contraction2 = {
          Eigen::IndexPair<Eigen::Index>(0, 0), Eigen::IndexPair<Eigen::Index>(1, 1)};
      const auto tCdual                  = tensorView(dualC, std::array<Eigen::Index, 2>({3, 3}));
      const auto tCdualT                 = tensorView(dualC2, std::array<Eigen::Index, 2>({3, 3}));
      const auto prod                    = Cmoduli.contract(tCdual, double_contraction);
      const Eigen::Tensor<double, 0> res = tCdualT.contract(prod, double_contraction2);
      e.grad.grad                        = res(0) / 4.0; // extracting value of zero order tensor

      return e;
    } else
      static_assert(Dune::AlwaysFalse<Derived>::value, "No fitting ScalarType.");
  }

  /**
   * \brief A helper function to compute the deviatoric part of the stresses.
   *
   * \details While using AutoDiff, if the eigenvalues (principal stretches) of C are degenerated, the derivative can
   * have certain singularities. This is circumvented here by explicitly updating the derivatives.
   *
   */
  template <typename Derived, typename ST>
  requires(std::same_as<typename Derived::Scalar, ST>)
  auto deviatoricStress(const Eigen::MatrixBase<Derived>& C, const Eigen::Vector<ST, dim>& lambdasST,
                        Eigen::Matrix<ST, dim, dim> NST) const {
    if constexpr (not Concepts::AutodiffScalar<ST>) {
      return transformDeviatoricStresses(dev_.stresses(lambdasST), NST);
    } else if constexpr (std::is_same_v<ST, autodiff::dual>) {
      constexpr int nVoigtIndices = 6;
      Eigen::Vector<autodiff::dual, nVoigtIndices> g;
      const auto Cvec         = toVoigt(C.derived());
      const auto realCVec     = derivative<0>(Cvec);
      const auto realC        = fromVoigt(realCVec);
      const auto dualC        = fromVoigt(Cvec.unaryExpr([](const auto& v) { return v.grad; }).eval());
      const auto [lambdas, N] = principalStretches(realC);

      const auto stresses = toVoigt(transformDeviatoricStresses(dev_.stresses(lambdas), N));
      const auto Cmoduli  = toVoigt(transformDeviatoricTangentModuli(dev_.tangentModuli(lambdas), N));
      Eigen::Vector<double, nVoigtIndices> stressDirectionalDerivatrive = Cmoduli * toVoigt(dualC);
      stressDirectionalDerivatrive.topRows<3>() /= 2.0;

      for (int i = 0; i < nVoigtIndices; ++i) {
        g[i].val  = stresses[i];
        g[i].grad = stressDirectionalDerivatrive[i];
      }

      return fromVoigt(g);
    } else
      static_assert(Dune::AlwaysFalse<Derived>::value, "No fitting ScalarType.");
  }
};

} // namespace Ikarus::Materials
