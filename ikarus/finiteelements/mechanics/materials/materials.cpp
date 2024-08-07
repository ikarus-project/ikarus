// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file materials.cpp
 * \brief Explicit instantiations for material templates.
 */

#include <ikarus/finiteelements/mechanics/materials/linearelasticity.hh>
#include <ikarus/finiteelements/mechanics/materials/neohooke.hh>
#include <ikarus/finiteelements/mechanics/materials/svk.hh>
#include <ikarus/finiteelements/mechanics/materials/vanishingstrain.hh>
#include <ikarus/finiteelements/mechanics/materials/vanishingstress.hh>

namespace Ikarus {

template struct LinearElasticityT<double>;
template struct NeoHookeT<double>;
template struct StVenantKirchhoffT<double>;

// Forward compile planeStress and planeStrain
template auto planeStress(const LinearElasticity& mat, typename LinearElasticity::ScalarType tol);
template auto planeStress(const NeoHooke& mat, typename NeoHooke::ScalarType tol);
template auto planeStress(const StVenantKirchhoff& mat, typename StVenantKirchhoff::ScalarType tol);

template auto planeStrain(const LinearElasticity& mat);
template auto planeStrain(const NeoHooke& mat);
template auto planeStrain(const StVenantKirchhoff& mat);

} // namespace Ikarus
