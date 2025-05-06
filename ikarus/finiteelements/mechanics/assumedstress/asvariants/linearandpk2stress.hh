// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file linearandpk2stress.hh
 * \brief Definition of the types of Assumed Stress formulations for 2D and 3D linear solid elements, where linear and
 * PK2 stresses are assumed.
 *
 * \ingroup ps
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::PS {

namespace Impl {
  template <typename GEO>
  auto transformationMatrixAtCenter(const GEO& geometry) {
    const auto& referenceElement = Dune::ReferenceElements<double, GEO::mydimension>::general(geometry.type());
    const auto quadPos0          = referenceElement.position(0, 0);

    return (transformationMatrix(geometry, quadPos0).transpose()).eval();
  }
} // namespace Impl

/**
 * \brief Interface for displacement-based Assumed Stress elements, where linear or PK2 stresses are assumed.
 *
 * \details See \cite viebahn_extension_2019 for details.
 */
template <typename GEO, int ass>
struct SX
{
  static constexpr int myDim             = GEO::mydimension;
  static constexpr int stressSize        = myDim * (myDim + 1) / 2;
  static constexpr int assumedStressSize = ass;
  using AnsatzType                       = Eigen::Matrix<double, stressSize, assumedStressSize>;
  using HType                            = Eigen::Matrix<double, assumedStressSize, assumedStressSize>;

  SX() = default;
  explicit SX(const GEO& geometry)
      : geometry_{std::make_optional<GEO>(geometry)},
        T0_{Impl::transformationMatrixAtCenter(geometry)} {}

protected:
  std::optional<GEO> geometry_;
  Eigen::Matrix<double, stressSize, stressSize> T0_;
};

/**
 * \brief S5 struct for AssumedStress (PS) for Q1 with 5 assumed stress parameters.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct S5 : SX<GEO, 5>
{
  using Base                             = SX<GEO, 5>;
  static constexpr int myDim             = Base::myDim;
  static constexpr int stressSize        = Base::stressSize;
  static constexpr int assumedStressSize = Base::assumedStressSize;
  using AnsatzType                       = typename Base::AnsatzType;
  using HType                            = typename Base::HType;

  S5() = default;
  explicit S5(const GEO& geo)
      : Base(geo) {}

  /**
   * \brief A function that implements the ansatz for the stress measure.
   * \details To align with conventions commonly found in the literature, where the integration domain is typically
   * [−1,1]^dim, a domain transformation is applied to convert from Dune's default domain of [0,1]^dim.
   */
  AnsatzType operator()(const Dune::FieldVector<double, myDim>& quadPos) const {
    AnsatzType P;
    P.setZero(stressSize, assumedStressSize);
    const double xi  = 2.0 * quadPos[0] - 1.0;
    const double eta = 2.0 * quadPos[1] - 1.0;
    P(0, 0)          = 1.0;
    P(1, 1)          = 1.0;
    P(2, 2)          = 1.0;
    P(0, 3)          = eta;
    P(1, 4)          = xi;

    return this->T0_ * P;
  }
};

/**
 * \brief S18 struct for AssumedStress (PS) for H1 with 18 stress parameters.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct S18 : SX<GEO, 18>
{
  using Base                             = SX<GEO, 18>;
  static constexpr int myDim             = Base::myDim;
  static constexpr int stressSize        = Base::stressSize;
  static constexpr int assumedStressSize = Base::assumedStressSize;
  using AnsatzType                       = typename Base::AnsatzType;
  using HType                            = typename Base::HType;

  S18() = default;
  explicit S18(const GEO& geo)
      : Base(geo) {}

  /**
   * \brief A function that implements the ansatz for the stress measure.
   * \details To align with conventions commonly found in the literature, where the integration domain is typically
   * [−1,1]^dim, a domain transformation is applied to convert from Dune's default domain of [0,1]^dim.
   */
  AnsatzType operator()(const Dune::FieldVector<double, myDim>& quadPos) const {
    AnsatzType P;
    P.setZero();
    const double xi   = 2.0 * quadPos[0] - 1.0;
    const double eta  = 2.0 * quadPos[1] - 1.0;
    const double zeta = 2.0 * quadPos[2] - 1.0;

    // ξξ = (1.0, eta, zeta, eta*zeta) → row 0:
    P(0, 0) = 1.0;
    P(0, 1) = eta;
    P(0, 2) = zeta;
    P(0, 3) = eta * zeta;

    // ηη = (1.0, xi, zeta, xi*zeta) → row 1:
    P(1, 4) = 1.0;
    P(1, 5) = xi;
    P(1, 6) = zeta;
    P(1, 7) = xi * zeta;

    // ζζ = (1.0, xi, eta, xi*eta) → row 2:
    P(2, 8)  = 1.0;
    P(2, 9)  = xi;
    P(2, 10) = eta;
    P(2, 11) = xi * eta;

    // ηζ = (1, xi) → row 3:
    P(3, 12) = 1.0;
    P(3, 13) = xi;

    // ξζ = (1.0, eta) → row 4:
    P(4, 14) = 1.0;
    P(4, 15) = eta;

    // ξη = (1.0, zeta) → row 5:
    P(5, 16) = 1.0;
    P(5, 17) = zeta;

    return this->T0_ * P;
  }
};

/**
 * \brief S24 struct for AssumedStress (PS) for H1 with 24 assumed stress parameters.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct S24 : SX<GEO, 24>
{
  using Base                             = SX<GEO, 24>;
  static constexpr int myDim             = Base::myDim;
  static constexpr int stressSize        = Base::stressSize;
  static constexpr int assumedStressSize = Base::assumedStressSize;
  using AnsatzType                       = typename Base::AnsatzType;
  using HType                            = typename Base::HType;

  S24() = default;
  explicit S24(const GEO& geo)
      : Base(geo) {}

  /**
   * \brief A function that implements the ansatz for the stress measure.
   * \details To align with conventions commonly found in the literature, where the integration domain is typically
   * [−1,1]^dim, a domain transformation is applied to convert from Dune's default domain of [0,1]^dim.
   */
  AnsatzType operator()(const Dune::FieldVector<double, myDim>& quadPos) const {
    AnsatzType P;
    P.setZero();
    const double xi   = 2.0 * quadPos[0] - 1.0;
    const double eta  = 2.0 * quadPos[1] - 1.0;
    const double zeta = 2.0 * quadPos[2] - 1.0;

    // ξξ = (1.0, eta, zeta, eta*zeta) → row 0:
    P(0, 0) = 1.0;
    P(0, 1) = eta;
    P(0, 2) = zeta;
    P(0, 3) = eta * zeta;

    // ηη = (1.0, xi, zeta, xi*zeta) → row 1:
    P(1, 4) = 1.0;
    P(1, 5) = xi;
    P(1, 6) = zeta;
    P(1, 7) = xi * zeta;

    // ζζ = (1.0, xi, eta, xi*eta) → row 2:
    P(2, 8)  = 1.0;
    P(2, 9)  = xi;
    P(2, 10) = eta;
    P(2, 11) = xi * eta;

    // ηζ = (1.0, xi, xi*eta, xi*zeta) → row 3:
    P(3, 12) = 1.0;
    P(3, 13) = xi;
    P(3, 14) = xi * eta;
    P(3, 15) = xi * zeta;

    // ξζ = (1.0, eta, xi*eta, eta*zeta) → row 4:
    P(4, 16) = 1.0;
    P(4, 17) = eta;
    P(4, 18) = xi * eta;
    P(4, 19) = eta * zeta;

    // ξη = (1.0, zeta, eta*zeta, xi*zeta) → row 5:
    P(5, 20) = 1.0;
    P(5, 21) = zeta;
    P(5, 22) = eta * zeta;
    P(5, 23) = xi * zeta;

    return this->T0_ * P;
  }
};

/**
 * \brief S30 struct for AssumedStress (PS) for H1 with 30 assumed stress parameters.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct S30 : SX<GEO, 30>
{
  using Base                             = SX<GEO, 30>;
  static constexpr int myDim             = Base::myDim;
  static constexpr int stressSize        = Base::stressSize;
  static constexpr int assumedStressSize = Base::assumedStressSize;
  using AnsatzType                       = typename Base::AnsatzType;
  using HType                            = typename Base::HType;

  S30() = default;
  explicit S30(const GEO& geo)
      : Base(geo) {}

  /**
   * \brief A function that implements the ansatz for the stress measure.
   * \details To align with conventions commonly found in the literature, where the integration domain is typically
   * [−1,1]^dim, a domain transformation is applied to convert from Dune's default domain of [0,1]^dim.
   */
  AnsatzType operator()(const Dune::FieldVector<double, myDim>& quadPos) const {
    AnsatzType P;
    P.setZero();
    const double xi   = 2.0 * quadPos[0] - 1.0;
    const double eta  = 2.0 * quadPos[1] - 1.0;
    const double zeta = 2.0 * quadPos[2] - 1.0;

    // ξξ = (1.0, eta, zeta, eta*zeta) → row 0:
    P(0, 0) = 1.0;
    P(0, 1) = eta;
    P(0, 2) = zeta;
    P(0, 3) = eta * zeta;

    // ηη = (1.0, xi, zeta, xi*zeta) → row 1:
    P(1, 4) = 1.0;
    P(1, 5) = xi;
    P(1, 6) = zeta;
    P(1, 7) = xi * zeta;

    // ζζ = (1.0, xi, eta, xi*eta) → row 2:
    P(2, 8)  = 1.0;
    P(2, 9)  = xi;
    P(2, 10) = eta;
    P(2, 11) = xi * eta;

    // ηζ = (1.0, xi, xi*eta, xi*zeta) → row 3:
    P(3, 12) = 1.0;
    P(3, 13) = xi;
    P(3, 14) = eta;
    P(3, 15) = zeta;
    P(3, 16) = xi * eta;
    P(3, 17) = xi * zeta;

    // ξζ = (1.0, eta, xi*eta, eta*zeta) → row 4:
    P(4, 18) = 1.0;
    P(4, 19) = xi;
    P(4, 20) = eta;
    P(4, 21) = zeta;
    P(4, 22) = xi * eta;
    P(4, 23) = eta * zeta;

    // ξη = (1.0, zeta, eta*zeta, xi*zeta) → row 5:
    P(5, 24) = 1.0;
    P(5, 25) = xi;
    P(5, 26) = eta;
    P(5, 27) = zeta;
    P(5, 28) = eta * zeta;
    P(5, 29) = xi * zeta;

    return this->T0_ * P;
  }
};

} // namespace Ikarus::PS
