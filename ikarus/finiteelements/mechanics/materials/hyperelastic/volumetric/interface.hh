// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Implementation of the volumetric part of a hyperelastic material.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Interface for the volumetric part of a hyperelastic material. Has to be parametrized with a volumetric
 * function.
 * \ingroup materials
 *
 * \details The volumetric part of the hyperelastic model, i.e., related to \f$ U(J) \f$, is
 * parametrized with a certain volumetric function (VF) implemented in terms of the determinant of the deformation
 * gradient. The underlying volumetric function must only implement the energy \f$ U(J) \f$ and its first
 * and second derivatives w.r.t \f$ J \f$.
 *
 * \tparam VF Volumetric function.
 */
template <Concepts::VolumetricFunction VF>
struct Volumetric
{
  using VolumetricFunction = VF;
  using MaterialParameter  = double;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Volumetric function: " + VF::name(); }

  /**
   * \brief Construct a new volumetric part
   *
   * \param matPar Materialparameter that is considered as the penalty parameter for the constraint against volumetric
   * deformations. Typically, if the deviatoric part of the energy function is a function of the total principal
   * stretches, then this material parameter is the  Lamé's first parameter and if the energy is a function of the
   * deviatoric principal stretches, then bulk modulus is used.
   *
   * \param vf the volumetric function
   */
  template <typename VFF>
  Volumetric(MaterialParameter matPar, VFF&& vf)
      : matPar_{matPar},
        volumetricFunction_{std::forward<VFF>(vf)} {}

  /**
   * \brief Returns the material parameters stored in the deviatoric part of the material.
   */
  const MaterialParameter materialParameter() const { return matPar_; }

  /**
   * \brief Computes stored energy of the volumetric function.
   *
   * \param J determinant of the deformation gradient \f$ J = \det\BF \f$.
   * \tparam ST the scalartype of J
   * \return ST energy.
   */
  template <typename ST>
  ST storedEnergy(const ST& J) const {
    return matPar_ * volumetricFunction_.storedEnergyImpl(J);
  };

  /**
   * \brief Computes the first derivatives of the energy of the volumetric function w.r.t \f$ J \f$.
   *
   * \param J determinant of the deformation gradient \f$ J = \det\BF \f$.
   * \tparam ST the scalartype of J
   * \return ST first derivative of the energy.
   */
  template <typename ST>
  ST firstDerivative(const ST& J) const {
    return matPar_ * volumetricFunction_.firstDerivativeImpl(J);
  }

  /**
   * \brief Computes the second derivatives of the energy of the volumetric function w.r.t \f$ J \f$.
   *
   * \param J determinant of the deformation gradient \f$ J = \det\BF \f$.
   * \tparam ST the scalartype of J
   * \return ST second derivative of the energy.
   */
  template <typename ST>
  ST secondDerivative(const ST& J) const {
    return matPar_ * volumetricFunction_.secondDerivativeImpl(J);
  };

private:
  MaterialParameter matPar_;
  VF volumetricFunction_;
};

#ifndef DOXYGEN
template <typename VF>
Volumetric(double, VF&&) -> Volumetric<std::remove_cvref_t<VF>>;
#endif

} // namespace Ikarus::Materials
