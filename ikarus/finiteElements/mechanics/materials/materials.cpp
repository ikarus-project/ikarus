
// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <ikarus/finiteElements/mechanics/materials/linearElasticity.hh>
#include <ikarus/finiteElements/mechanics/materials/neohooke.hh>
#include <ikarus/finiteElements/mechanics/materials/svk.hh>
namespace Ikarus {

  template struct LinearElasticityT<double>;
  template struct NeoHookeT<double>;
  template struct StVenantKirchhoffT<double>;
}  // namespace Ikarus
