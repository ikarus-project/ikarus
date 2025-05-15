// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcantileverbeam.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/assumedstress.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mixin.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <typename MAT>
auto cantileverBeamResults(const MAT& mat) {
  static_assert(MAT::isReduced,
                "cantileverBeamResults are available only for a reduced material (planeStress or planeStrain).");
  using namespace Materials;
  using MATU = typename MAT::Underlying;

  if constexpr (std::is_same_v<MATU, StVenantKirchhoff>)
    return std::make_pair(60, 4.461888037704155);
  else if constexpr (std::is_same_v<MATU, NeoHooke>)
    return std::make_pair(60, 4.481945889821968);
  else if constexpr (std::is_same_v<MATU, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>)
    return std::make_pair(60, 4.6903421642761485);
  else
    static_assert(Dune::AlwaysFalse<MATU>::value, "Expected results are not available for the given material.");
}

template <typename MAT>
auto cantileverBeamResults3D(const MAT& mat, int numberOfParameters) {
  using namespace Materials;

  if (numberOfParameters == 18) {
    if constexpr (std::is_same_v<MAT, StVenantKirchhoff>)
      return std::make_pair(60, 4.772076885669366);
    else if constexpr (std::is_same_v<MAT, NeoHooke>)
      return std::make_pair(60, 4.784585421202026);
    else if constexpr (std::is_same_v<MAT, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>)
      return std::make_pair(60, 4.9200507102060955);
    else
      static_assert(Dune::AlwaysFalse<MAT>::value, "Expected results are not available for the given material.");
  } else if (numberOfParameters == 24) {
    if constexpr (std::is_same_v<MAT, StVenantKirchhoff>)
      return std::make_pair(60, 4.751832972951744);
    else if constexpr (std::is_same_v<MAT, NeoHooke>)
      return std::make_pair(60, 4.765489302035456);
    else if constexpr (std::is_same_v<MAT, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>)
      return std::make_pair(60, 4.906056299966376);
    else
      static_assert(Dune::AlwaysFalse<MAT>::value, "Expected results are not available for the given material.");
  } else if (numberOfParameters == 30) {
    if constexpr (std::is_same_v<MAT, StVenantKirchhoff>)
      return std::make_pair(60, 4.370359429120268);
    else if constexpr (std::is_same_v<MAT, NeoHooke>)
      return std::make_pair(60, 4.376829591567893);
    else if constexpr (std::is_same_v<MAT, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>)
      return std::make_pair(60, 4.500152745952108);
    else
      static_assert(Dune::AlwaysFalse<MAT>::value, "Expected results are not available for the given material.");
  } else
    assert("Expected results are not available for the given number of parameters");
  __builtin_unreachable();
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameter);
  Materials::NeoHooke matNH(matParameter);
  auto matBK = Materials::makeBlatzKo(40.0);

  auto reducedMats = Dune::makeTupleVector(planeStrain(matSVK), planeStrain(matNH), planeStrain(matBK));
  Dune::Hybrid::forEach(reducedMats, [&t](const auto& mat) {
    t.subTest(
        cantileverBeamTest<2>(mat, skills(assumedStress<PS::PK2Stress>(5)), cantileverBeamResults(mat), false, true));
  });

  auto materials = Dune::makeTupleVector(matSVK, matNH, matBK);
  Dune::Hybrid::forEach(materials, [&t](const auto& mat) {
    for (int parameters : std::array{18, 24, 30})
      t.subTest(cantileverBeamTest<3>(mat, skills(assumedStress<PS::PK2Stress>(parameters)),
                                      cantileverBeamResults3D(mat, parameters)));
  });

  return t.exit();
}
