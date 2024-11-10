// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <muesli/muesli.h>

#include <ikarus/finiteelements/physicshelper.hh>

namespace Ikarus::Materials::Muesli {

// Alias for Muesli materials

using LinearElasticity = muesli::elasticAnisotropicMaterial;
using SVK              = muesli::svkMaterial;
using NeoHooke         = muesli::neohookeanMaterial;
using MooneyRivlin     = muesli::mooneyMaterial;
using Yeoh             = muesli::yeohMaterial;
using ArrudaBoyce      = muesli::arrudaboyceMaterial;

// Alias for Muesli material properties
using MaterialProperties = muesli::materialProperties;

inline MaterialProperties propertiesFromIkarusMaterialParameters(const LamesFirstParameterAndShearModulus& mpt,
                                                                 bool regularized = false) {
  auto mpm = muesli::materialProperties{};
  mpm.insert({"lambda", mpt.lambda});
  mpm.insert({"mu", mpt.mu});
  if (regularized)
    mpm.insert({"subtype regularized", 0});

  return mpm;
}

inline MaterialProperties propertiesFromIkarusMaterialParameters(const YoungsModulusAndPoissonsRatio& mpt,
                                                                 bool regularized = false) {
  auto mpm = muesli::materialProperties{};
  mpm.insert({"young", mpt.emodul});
  mpm.insert({"poisson", mpt.nu});

  if (regularized)
    mpm.insert({"subtype regularized", 0});

  return mpm;
}

} // namespace Ikarus::Materials::Muesli