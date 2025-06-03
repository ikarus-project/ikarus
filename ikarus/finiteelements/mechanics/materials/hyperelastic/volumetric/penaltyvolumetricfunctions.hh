// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file volumetricfunctions.hh
 * \brief Implementation of the volumetric part of a hyperelastic material.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/volumetric/interface.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

struct PVF1
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return J - 1;
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& /* J */) const {
    return 1;
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& /* J */) const {
    return 0.0;
  };

  [[nodiscard]] constexpr static std::string name() noexcept { return "Penalty Function 1"; }
};

} // namespace Ikarus::Materials
