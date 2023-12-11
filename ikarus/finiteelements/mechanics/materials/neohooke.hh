// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

  // see Javier Bonet, Richard D. Wood - Nonlinear Continuum Mechanics for Finite Element Analysis, 2nd Edition (2008)
  // Section 6.4.3
  template <typename ScalarType_ = double>
  struct NeoHookeT : public Material<NeoHookeT<ScalarType_>> {
    [[nodiscard]] constexpr std::string nameImpl() const noexcept { return "NeoHooke"; }

    explicit NeoHookeT(const LamesFirstParameterAndShearModulus& mpt) : lambdaAndmu{mpt} {}
    using ScalarType                    = ScalarType_;
    static constexpr int worldDimension = 3;
    using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
    using StressMatrix                  = StrainMatrix;

    static constexpr auto strainTag          = StrainTags::rightCauchyGreenTensor;
    static constexpr auto stressTag          = StressTags::PK2;
    static constexpr auto tangentModuliTag   = TangentModuliTags::Material;
    static constexpr bool energyAcceptsVoigt = false;
    static constexpr bool stressToVoigt      = false;
    static constexpr bool stressAcceptsVoigt = false;
    static constexpr bool moduliToVoigt      = false;
    static constexpr bool moduliAcceptsVoigt = false;
    // this factor denotes the differences between the returned stresses and moduli and the passed strain
    // for neoHooke the inserted quantity is C the right Cauchy-Green tensor,
    // the function relation between the energy and the stresses is S = 2*\partial \psi(C)/ \partial C.
    // This factor is pre factor, which is the difference to the actual derivative is written here
    static constexpr double derivativeFactor = 2;

    template <typename Derived>
    ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& C) const noexcept {
      static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
      if constexpr (!Concepts::EigenVector<Derived>) {
        const auto traceC  = C.trace();
        const auto logdetF = log(sqrt(C.determinant()));
        return lambdaAndmu.mu / 2.0 * (traceC - 3 - 2 * logdetF) + lambdaAndmu.lambda / 2.0 * logdetF * logdetF;
      } else
        static_assert(!Concepts::EigenVector<Derived>,
                      "NeoHooke energy can only be called with matrix and not a vector in Voigt notation");
    }

    template <bool voigt, typename Derived>
    auto stressesImpl(const Eigen::MatrixBase<Derived>& C) const {
      static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
      if constexpr (!voigt) {
        if constexpr (!Concepts::EigenVector<Derived>) {
          const auto logdetF = log(sqrt(C.determinant()));
          const auto invC    = C.inverse().eval();
          return (lambdaAndmu.mu * (StrainMatrix::Identity() - invC) + lambdaAndmu.lambda * logdetF * invC).eval();
        } else
          static_assert(!Concepts::EigenVector<Derived>,
                        "NeoHooke can only be called with matrix and not a vector in Voigt notation");
      } else
        static_assert(voigt == false, "NeoHooke does not support returning stresses in Voigt notation");
    }

    template <bool voigt, typename Derived>
    auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& C) const {
      static_assert(Concepts::EigenMatrixOrVoigtNotation3<Derived>);
      if constexpr (!voigt) {
        const auto invC    = C.inverse().eval();
        const auto logdetF = log(sqrt(C.determinant()));
        const auto CTinv   = TensorCast(invC, std::array<Eigen::Index, 2>({3, 3}));
        static_assert(Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3>>::NumIndices == 2);
        // see Javier Bonet, Richard D. Wood - Nonlinear Continuum Mechanics for Finite Element Analysis, 2nd Edition
        // (2008) Eq. 6.30
        //  (first edition has an error)
        Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<3, 3, 3, 3>> moduli
            = (lambdaAndmu.lambda * dyadic(CTinv, CTinv)
               + 2 * (lambdaAndmu.mu - lambdaAndmu.lambda * logdetF) * symTwoSlots(fourthOrderIKJL(invC, invC), {2, 3}))
                  .eval();
        return moduli;
      } else
        static_assert(voigt == false, "NeoHooke does not support returning tangentModuli in Voigt notation");
    }

    template <typename ScalarTypeOther>
    auto rebind() const {
      return NeoHookeT<ScalarTypeOther>(lambdaAndmu);
    }

    LamesFirstParameterAndShearModulus lambdaAndmu;
  };
  typedef NeoHookeT<double> NeoHooke;

}  // namespace Ikarus
