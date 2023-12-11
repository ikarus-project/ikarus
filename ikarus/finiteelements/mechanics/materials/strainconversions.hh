// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once
#include "tags.hh"

#include <unsupported/Eigen/MatrixFunctions>

#include <Eigen/Core>

#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {
  namespace Impl {
    template <typename Derived>
    decltype(auto) mayTransformFromVoigt(const Eigen::MatrixBase<Derived>& e, bool isStrain = true) {
      if constexpr (Concepts::EigenVector<Derived>)
        return fromVoigt(e.derived(), isStrain);
      else
        return e;
    }
  }  // namespace Impl

  template <StrainTags tag, typename Derived>
  auto createGreenLagrangianStrains(const Eigen::MatrixBase<Derived>& eMB) {
    const auto& e = eMB.derived();
    static_assert(Concepts::EigenMatrix33<Derived>);
    if constexpr (tag == StrainTags::greenLagrangian)
      return e;
    else if constexpr (tag == StrainTags::deformationGradient)
      return (0.5 * (e.transpose() * e - Derived::Identity())).eval();
    else if constexpr (tag == StrainTags::displacementGradient)
      return (0.5 * (e + e.transpose() + e.transpose() * e)).eval();
    else if constexpr (tag == StrainTags::rightCauchyGreenTensor)
      return (0.5 * (e - Derived::Identity())).eval();
  }

  template <StrainTags tag, typename Derived>
  decltype(auto) createDeformationGradient(const Eigen::MatrixBase<Derived>& eMB) {
    const auto& e = eMB.derived();

    static_assert(Concepts::EigenMatrix33<Derived>);
    if constexpr (tag == StrainTags::greenLagrangian) {
      // E = 0.5 * (F ^ 2 - I);
      // 2*E = F ^ 2 - I;
      // 2*E+I = F ^ 2;
      // sqrt(2*E+I) = F;
      return ((2 * e + Derived::Identity()).sqrt()).eval();
    } else if constexpr (tag == StrainTags::deformationGradient)
      return e;
    else if constexpr (tag == StrainTags::displacementGradient)
      return (e + Derived::Identity()).eval();
    else if constexpr (tag == StrainTags::rightCauchyGreenTensor) {
      return (e.sqrt()).eval();  // this looses information, since the rotation information from an original F is lost!
    }
  }

  template <StrainTags tag, typename Derived>
  decltype(auto) createRightCauchyGreen(const Eigen::MatrixBase<Derived>& eMB) {
    const auto& e = eMB.derived();
    static_assert(Concepts::EigenMatrix33<Derived>);
    if constexpr (tag == StrainTags::greenLagrangian) {
      // E = 0.5 * (C - I);
      // 2*E = C - I;
      // 2*E+I = C;
      return (2 * e + Derived::Identity()).eval();
    } else if constexpr (tag == StrainTags::deformationGradient)
      return (e.transpose() * e).eval();
    else if constexpr (tag == StrainTags::displacementGradient) {
      const auto F = e + Derived::Identity();
      return (F.transpose() * F).eval();
    } else if constexpr (tag == StrainTags::rightCauchyGreenTensor) {
      return e;
    }
  }

  template <StrainTags from, StrainTags to, typename Derived>
  decltype(auto) transformStrain(const Eigen::MatrixBase<Derived>& eRaw) {
    static_assert((from == to) or (from != StrainTags::linear and to != StrainTags::linear),
                  "No useful transformation available for linear strains.");
    decltype(auto) e = Impl::mayTransformFromVoigt(eRaw, true);
    if constexpr (from == to)
      return e;
    else if constexpr (to == StrainTags::greenLagrangian)
      return createGreenLagrangianStrains<from>(e);
    else if constexpr (to == StrainTags::deformationGradient)
      return createDeformationGradient<from>(e);
    else if constexpr (to == StrainTags::rightCauchyGreenTensor) {
      return createRightCauchyGreen<from>(e);
    } else
      static_assert(to == StrainTags::greenLagrangian or to == StrainTags::deformationGradient
                    or to == StrainTags::rightCauchyGreenTensor);
  }
}  // namespace Ikarus
