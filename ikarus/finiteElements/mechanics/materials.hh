// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteElements/mechanics/materials/interface.hh>
#include <ikarus/finiteElements/mechanics/materials/linearElasticity.hh>
#include <ikarus/finiteElements/mechanics/materials/neohooke.hh>
#include <ikarus/finiteElements/mechanics/materials/svk.hh>
#include <ikarus/finiteElements/mechanics/materials/tags.hh>
#include <ikarus/finiteElements/mechanics/materials/vanishingStress.hh>

namespace Ikarus {
  template <typename MaterialImpl>
  auto planeStress(const MaterialImpl& mat, typename MaterialImpl::ScalarType p_tol = 1e-8) {
    return makeVanishingStress<{2, 1}, {2, 0}, {2, 2}>(mat, p_tol);
  }

  template <typename MaterialImpl>
  auto shellMaterial(const MaterialImpl& mat, typename MaterialImpl::ScalarType p_tol = 1e-8) {
    return makeVanishingStress<{2, 2}>(mat, p_tol);
  }

  template <typename MaterialImpl>
  auto beamMaterial(const MaterialImpl& mat, typename MaterialImpl::ScalarType p_tol = 1e-8) {
    return makeVanishingStress<{1, 1}, {2, 2}>(mat, p_tol);
  }
}  // namespace Ikarus
