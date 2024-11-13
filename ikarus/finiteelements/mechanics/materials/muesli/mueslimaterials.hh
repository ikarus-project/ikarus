// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.hh
 * \brief Header file for material models in Ikarus finite element mechanics.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/muesli/muesliefinite.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslielastic.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslihelpers.hh>

namespace Ikarus::Materials::Muesli {

// Alias for Muesli materials

using LinearElasticity            = muesli::elasticIsotropicMaterial;
using LinearAnisotropicElasticity = muesli::elasticAnisotropicMaterial;
using LinearOrthotropicElasticity = muesli::elasticOrthotropicMaterial;

using StVenantKirchhoff = muesli::svkMaterial;
using NeoHooke          = muesli::neohookeanMaterial;
using MooneyRivlin      = muesli::mooneyMaterial;
using Yeoh              = muesli::yeohMaterial;
using ArrudaBoyce       = muesli::arrudaboyceMaterial;

template <typename MPT>
requires(std::same_as<MPT, YoungsModulusAndPoissonsRatio> or std::same_as<MPT, LamesFirstParameterAndShearModulus>)
auto makeNeoHooke(const MPT& mpt, bool useDeviatoricStretches = false) {
  auto muesliParameters = propertiesFromIkarusMaterialParameters(mpt);
  if (useDeviatoricStretches)
    addRegularizedTag(muesliParameters);
  return MuesliFinite<Muesli::NeoHooke>(muesliParameters);
}

template <typename MPT>
requires(std::same_as<MPT, YoungsModulusAndPoissonsRatio> or std::same_as<MPT, LamesFirstParameterAndShearModulus>)
auto makeSVK(const MPT& mpt) {
  auto muesliParameters = propertiesFromIkarusMaterialParameters(mpt);
  return MuesliFinite<Muesli::StVenantKirchhoff>(muesliParameters);
}

} // namespace Ikarus::Materials::Muesli