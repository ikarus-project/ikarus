// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file displacementgradient.hh
 * \brief Definition of the types of EAS formulations for 2D and 3D linear solid elements, where the displacement
 * gradient is enhanced.
 *
 * \ingroup eas
 */

#pragma once

#include <ikarus/finiteelements/mechanics/strainenhancements/easvariants/helperfunctions.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS {

/**
 * \brief Interface for displacement-based EAS elements, where displacement gradient is enhanced.
 *
 * \details See \cite pfefferkornTransformationsShapeFunctions2019b for details.
 */
template <typename GEO, int ess>
struct HX
{
  static constexpr int myDim              = GEO::mydimension;
  static constexpr int strainSize         = myDim * (myDim + 1) / 2;
  static constexpr int enhancedStrainSize = ess;
  using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;
  using DType                             = Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>;
  using AnsatzType                        = Eigen::Matrix<double, myDim, myDim>;
  using HType                             = std::array<AnsatzType, enhancedStrainSize>;

  HX() = default;
  explicit HX(const GEO& geometry)
      : geometry_{std::make_optional<GEO>(geometry)} {}

protected:
  std::optional<GEO> geometry_;
  const auto& geometry() const { return this->geometry_.value(); };
};

/**
 * \brief Dummy struct for displacement-based EAS elements, i.e. 0 enhanced modes, where displacement gradient is
 * enhanced.
 */
template <typename GEO>
struct H0 : HX<GEO, 0>
{
  using Base                              = HX<GEO, 0>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;
  using AnsatzType                        = typename Base::AnsatzType;
  using HType                             = typename Base::HType;

  H0() = default;
  explicit H0(const GEO& geo)
      : Base(geo) {}

  HType operator()(const Dune::FieldVector<double, myDim>& /*quadPos*/) const { return {}; }
};

/**
 * \brief H4 struct for EAS with 4 enhanced modes.
 *
 * \details The H4 struct represents an implementation of EAS for a specific case where displacement gradient is
 * enhanced (2D case).
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct H4 : HX<GEO, 4>
{
  using Base                              = HX<GEO, 4>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;
  using AnsatzType                        = typename Base::AnsatzType;
  using HType                             = typename Base::HType;

  H4() = default;
  explicit H4(const GEO& geo)
      : Base(geo) {}

  HType operator()(const Dune::FieldVector<double, myDim>& quadPos) const {
    const double xi    = quadPos[0];
    const double eta   = quadPos[1];
    const double xi_t  = 2 * xi - 1.0;
    const double eta_t = 2 * eta - 1.0;

    HType H;
    std::ranges::fill(H, AnsatzType::Zero());

    H[0](0, 0) = xi_t;
    H[1](0, 1) = eta_t;
    H[2](1, 0) = xi_t;
    H[3](1, 1) = eta_t;

    for (auto& HMat : H)
      Impl::transformDisplacementGradient(this->geometry(), HMat, quadPos);

    return H;
  }
};

/**
 * \brief H9 struct for EAS with 9 enhanced modes.
 *
 * \details The H9 struct represents an implementation of EAS for a specific case where displacement gradient is
 * enhanced (3D case).
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct H9 : HX<GEO, 9>
{
  using Base                              = HX<GEO, 9>;
  static constexpr int myDim              = Base::myDim;
  static constexpr int strainSize         = Base::strainSize;
  static constexpr int enhancedStrainSize = Base::enhancedStrainSize;
  using MType                             = typename Base::MType;
  using DType                             = typename Base::DType;
  using AnsatzType                        = typename Base::AnsatzType;
  using HType                             = typename Base::HType;

  H9() = default;
  explicit H9(const GEO& geo)
      : Base(geo) {}

  HType operator()(const Dune::FieldVector<double, myDim>& quadPos) const {
    const double xi     = quadPos[0];
    const double eta    = quadPos[1];
    const double zeta   = quadPos[2];
    const double xi_t   = 2 * xi - 1.0;
    const double eta_t  = 2 * eta - 1.0;
    const double zeta_t = 2 * zeta - 1.0;

    HType H;
    std::ranges::fill(H, AnsatzType::Zero());

    H[0](0, 0) = xi_t;
    H[1](0, 1) = eta_t;
    H[2](0, 2) = zeta_t;
    H[3](1, 0) = xi_t;
    H[4](1, 1) = eta_t;
    H[5](1, 2) = zeta_t;
    H[6](2, 0) = xi_t;
    H[7](2, 1) = eta_t;
    H[8](2, 2) = zeta_t;

    for (auto& HMat : H)
      Impl::transformDisplacementGradient(this->geometry(), HMat, quadPos);

    return H;
  }
};
} // namespace Ikarus::EAS
