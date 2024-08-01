// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ikarus/utils/makeenum.hh>

#pragma once
namespace Ikarus::Vtk {
/**
 * \enum DataTag
 * \brief Tag structure indicating cell data or point data.
 */
MAKE_ENUM(DataTag, asCellData, asPointData);
} // namespace Ikarus::Vtk