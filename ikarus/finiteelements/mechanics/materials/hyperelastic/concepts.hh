// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file concepts.hh
 * \brief Header file including concepts for hyperelastic material models.
 * \ingroup  materials
 */

#pragma once

#include <concepts>
#include <string>

#include <ikarus/utils/concepts.hh>

namespace Ikarus::Concepts {

/**
 * \concept DeviatoricFunction
 * \brief Concept to check if the underlying function is a deviatoric function.
 * \tparam Type of the deviatoric function.
 *
 */
template <typename DF>
concept DeviatoricFunction = requires(DF dm, const typename DF::template PrincipalStretches<>& lambda) {
  typename DF::ScalarType;

  typename DF::template PrincipalStretches<>;
  typename DF::template FirstDerivative<>;
  typename DF::template SecondDerivative<>;

  typename DF::MaterialParameters;

  requires Concepts::EigenVector<typename DF::template FirstDerivative<>>;
  requires Concepts::EigenMatrix<typename DF::template SecondDerivative<>>;
  requires std::is_same_v<typename DF::template PrincipalStretches<>, typename DF::template FirstDerivative<>>;

  { dm.storedEnergyImpl(lambda) } -> std::same_as<typename DF::ScalarType>;
  { dm.firstDerivativeImpl(lambda) } -> std::same_as<typename DF::template FirstDerivative<>>;
  { dm.secondDerivativeImpl(lambda) } -> std::same_as<typename DF::template SecondDerivative<>>;
  { dm.materialParametersImpl() } -> std::same_as<typename DF::MaterialParameters>;
  { dm.name() } -> std::convertible_to<std::string>;
};

/**
 * \concept VolumetricFunction
 * \brief Concept to check if the underlying function is a volumetric function.
 * \tparam Type of the volumetric function.
 *
 */
template <typename VF>
concept VolumetricFunction = requires(VF vf, const double& j) {
  { vf.storedEnergyImpl(j) } -> std::same_as<double>;
  { vf.firstDerivativeImpl(j) } -> std::same_as<double>;
  { vf.secondDerivativeImpl(j) } -> std::same_as<double>;
  { vf.name() } -> std::convertible_to<std::string>;
};

} // namespace Ikarus::Concepts
