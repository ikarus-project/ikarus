// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/utils/refl.hpp>

namespace Ikarus {

struct FESettingsContainer
{
  int nGP            = std::numeric_limits<int>::max();
  int orderGP        = std::numeric_limits<int>::max();
  float someSettings = std::numeric_limits<float>::max();
};

} // namespace Ikarus

// clang-format off
REFL_AUTO(
  type(Ikarus::FESettingsContainer),
  field(nGP),
  field(orderGP)
  )

