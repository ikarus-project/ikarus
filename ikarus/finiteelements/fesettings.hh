// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/utils/refl.hpp>

namespace Ikarus {

// Macro definition
#define registerSetting(name, type) type name{};

struct SettingsContainer
{
  registerSetting(nGP, int);
  registerSetting(orderGP, int);
};

} // namespace Ikarus

REFL_AUTO(type(Ikarus::SettingsContainer), field(nGP), field(orderGP))