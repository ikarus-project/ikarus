// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <muesli/muesli.h>

#include <ikarus/finiteelements/physicshelper.hh>

namespace Ikarus::Materials::Muesli {

// Alias for Muesli material properties
using MaterialProperties = muesli::materialProperties;

inline MaterialProperties propertiesFromIkarusMaterialParameters(const LamesFirstParameterAndShearModulus& mpt) {
  auto mpm = muesli::materialProperties{};
  mpm.insert({"lambda", mpt.lambda});
  mpm.insert({"mu", mpt.mu});

  return mpm;
}

inline void addRegularizedTag(MaterialProperties& mpm) { mpm.insert({"subtype regularized", 0}); }

inline MaterialProperties propertiesFromIkarusMaterialParameters(const YoungsModulusAndPoissonsRatio& mpt) {
  auto mpm = muesli::materialProperties{};
  mpm.insert({"young", mpt.emodul});
  mpm.insert({"poisson", mpt.nu});

  return mpm;
}

} // namespace Ikarus::Materials::Muesli