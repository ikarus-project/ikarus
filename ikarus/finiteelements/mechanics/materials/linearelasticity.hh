// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "svk.hh"

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/utils/tensorutils.hh>
namespace Ikarus {
  template <typename ScalarType_ = double>
  struct LinearElasticityT : public Material<LinearElasticityT<ScalarType_>> {
    [[nodiscard]] constexpr std::string nameImpl() const noexcept { return "LinearElasticity"; }

    using ScalarType = ScalarType_;
    using Base       = StVenantKirchhoffT<ScalarType>;

    explicit LinearElasticityT(const LamesFirstParameterAndShearModulus& mpt) : svk{mpt} {}
    using field_type                    = ScalarType;
    static constexpr int worldDimension = 3;
    using StrainMatrix                  = Eigen::Matrix<ScalarType, worldDimension, worldDimension>;
    using StressMatrix                  = StrainMatrix;

    static constexpr auto strainTag          = StrainTags::linear;
    static constexpr auto stressTag          = StressTags::linear;
    static constexpr auto tangentModuliTag   = TangentModuliTags::Material;
    static constexpr bool energyAcceptsVoigt = Base::energyAcceptsVoigt;
    static constexpr bool stressToVoigt      = Base::stressToVoigt;
    static constexpr bool stressAcceptsVoigt = Base::stressAcceptsVoigt;
    static constexpr bool moduliToVoigt      = Base::moduliToVoigt;
    static constexpr bool moduliAcceptsVoigt = Base::moduliAcceptsVoigt;
    // this factor denotes the differences between the returned stresses and moduli and the passed strain
    // for neoHooke the inserted quantity is C the right Cauchy-Green tensor,
    // the function relation between the energy and the stresses is S = 2*\partial \psi(C)/ \partial C.
    // This factor is pre factor, which is the difference to the actual derivative is written here
    static constexpr double derivativeFactor = 1;

    template <typename Derived>
    ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived>& E) const {
      return svk.template storedEnergyImpl(E);
    }

    template <bool voigt, typename Derived>
    auto stressesImpl(const Eigen::MatrixBase<Derived>& E) const {
      return svk.template stressesImpl<voigt>(E);
    }

    template <bool voigt, typename Derived>
    auto tangentModuliImpl(const Eigen::MatrixBase<Derived>& E) const {
      return svk.template tangentModuliImpl<voigt>(E);
    }

    template <typename ScalarTypeOther>
    auto rebind() const {
      LinearElasticityT<ScalarTypeOther> leRebound;
      leRebound.svk = svk.template rebind<ScalarTypeOther>();
      return leRebound;
    }

    StVenantKirchhoffT<ScalarType> svk;
  };

  typedef LinearElasticityT<double> LinearElasticity;
}  // namespace Ikarus
