// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {
  template <typename ScalarType_ = double>
  struct StVenantKirchhoffT : public Material<StVenantKirchhoffT<ScalarType_>> {
    [[nodiscard]] constexpr std::string nameImpl() const { return "StVenantKirchhoff"; }

    explicit StVenantKirchhoffT(const LamesFirstParameterAndShearModulus& mpt) : materialParameter{mpt} {}
    using ScalarType                    = ScalarType_;
    static constexpr int worldDimension = 3;
    using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
    using StressMatrix                  = StrainMatrix;

    static constexpr auto strainTag          = StrainTags::greenLagrangian;
    static constexpr auto stressTag          = StressTags::PK2;
    static constexpr auto tangentModuliTag   = TangentModuliTags::Material;
    static constexpr bool energyAcceptsVoigt = true;
    static constexpr bool stressToVoigt      = true;
    static constexpr bool stressAcceptsVoigt = true;
    static constexpr bool moduliToVoigt      = true;
    static constexpr bool moduliAcceptsVoigt = true;
    // this factor denotes the differences between the returned stresses and moduli and the passed strain
    // for neoHooke the inserted quantity is C the Green-Lagrangian strain tensor,
    // the function relation between the energy and the stresses is S = 1\partial \psi(E)/ \partial E.
    // This factor is pre factor, which is the difference to the actual derivative is written here
    static constexpr double derivativeFactor = 1;

    template <typename Derived>
    ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
      static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
      if constexpr (Concepts::EigenVector<Derived>) {
        const ScalarType traceE = E.template segment<3>(0).sum();
        const ScalarType squaredNorm
            = E.template segment<3>(0).squaredNorm() + E.template segment<3>(3).squaredNorm() / ScalarType(2.0);
        return materialParameter.lambda / ScalarType(2.0) * traceE * traceE + materialParameter.mu * squaredNorm;
      } else {
        const auto traceE = E.trace();
        return materialParameter.lambda / ScalarType(2.0) * traceE * traceE + materialParameter.mu * E.squaredNorm();
      }
    }

    template <bool voigt, typename Derived>
    auto stressesImpl(const Eigen::MatrixBase<Derived>& EMB) const {
      static_assert(Concepts::EigenMatrixOrVoigtNotation3<decltype(EMB.eval())>);
      const auto& E = EMB.derived();
      if constexpr (!voigt) {
        if constexpr (Concepts::EigenVector<Derived>) {
          static_assert(Concepts::EigenVector6<Derived>);
          Eigen::Matrix<ScalarType, 3, 3> S;
          const ScalarType traceE = E.template segment<3>(0).sum();
          S.diagonal().array()
              = materialParameter.lambda * traceE + 2 * materialParameter.mu * E.template segment<3>(0).array();
          S(1, 2) = S(2, 1) = materialParameter.mu * E(3);
          S(0, 2) = S(2, 0) = materialParameter.mu * E(4);
          S(0, 1) = S(1, 0) = materialParameter.mu * E(5);  // no two since E voigt has 2* on off-diagonal terms
          return S;
        } else {
          static_assert(Concepts::EigenMatrix33<Derived>);
          return (materialParameter.lambda * E.trace() * StrainMatrix::Identity() + 2 * materialParameter.mu * E)
              .eval();
        }
      } else {
        if constexpr (Concepts::EigenVector<Derived>) {
          static_assert(Concepts::EigenVector6<Derived>);
          Eigen::Matrix<ScalarType, 6, 1> S;
          const ScalarType traceE          = E.template segment<3>(0).sum();
          S.template segment<3>(0).array() = traceE * materialParameter.lambda;
          S.template segment<3>(0) += materialParameter.mu * 2 * E.template segment<3>(0);
          S.template segment<3>(3)
              = materialParameter.mu * E.template segment<3>(3);  // no two since E voigt has 2* on off-diagonal terms
          return S;
        } else {
          Eigen::Matrix<ScalarType, 6, 1> S;
          S.template segment<3>(0).array() = E.trace() * materialParameter.lambda;
          S.template segment<3>(0) += 2 * materialParameter.mu * E.diagonal();
          S(3) = 2 * materialParameter.mu * E(1, 2);
          S(4) = 2 * materialParameter.mu * E(0, 2);
          S(5) = 2 * materialParameter.mu * E(0, 1);
          return S;
        }
      }
    }

    template <bool voigt, typename Derived>
    auto tangentModuliImpl(const Eigen::MatrixBase<Derived>&) const {
      static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
      if constexpr (!voigt) {
        Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> moduli;
        moduli = materialParameter.lambda * identityFourthOrder()
                 + 2 * materialParameter.mu * symmetricIdentityFourthOrder();
        return moduli;
      } else {
        Eigen::Matrix<ScalarType, 6, 6> moduli;
        moduli.setZero();
        moduli.template block<3, 3>(0, 0).array() = materialParameter.lambda;
        moduli.template block<3, 3>(0, 0).diagonal().array() += 2 * materialParameter.mu;
        moduli.template block<3, 3>(3, 3).diagonal().array() = materialParameter.mu;
        return moduli;
      }
    }

    template <typename ScalarTypeOther>
    auto rebind() const {
      return StVenantKirchhoff<ScalarTypeOther>(materialParameter);
    }

    LamesFirstParameterAndShearModulus materialParameter;
  };

  typedef StVenantKirchhoffT<double> StVenantKirchhoff;

}  // namespace Ikarus
