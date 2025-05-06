// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcantileverbeam.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mixin.hh>
using namespace Ikarus;
using Dune::TestSuite;

template <typename ES, typename MAT>
auto cantileverBeamResults(const MAT& mat) {
  static_assert(MAT::isReduced,
                "cantileverBeamResults are available only for a reduced material (planeStress or planeStrain).");
  using namespace Materials;
  using MATU = typename MAT::Underlying;
  if constexpr (std::same_as<ES, EAS::GreenLagrangeStrain>) {
    if constexpr (std::is_same_v<MATU, StVenantKirchhoff>) {
      return std::make_pair(80, 4.459851990257645);
    } else if constexpr (std::is_same_v<MATU, NeoHooke>) {
      return std::make_pair(80, 4.479930218997457);
    } else if constexpr (std::is_same_v<MATU, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>) {
      return std::make_pair(80, 4.6877262928164365);
    } else
      static_assert(Dune::AlwaysFalse<MATU>::value,
                    "Expected results are not available for the given material and enhancedStrain type.");
  } else if constexpr (std::same_as<ES, EAS::DisplacementGradient>) {
    if constexpr (std::is_same_v<MATU, StVenantKirchhoff>) {
      return std::make_pair(80, 4.466045705106845);
    } else if constexpr (std::is_same_v<MATU, NeoHooke>) {
      return std::make_pair(80, 4.484284943716379);
    } else if constexpr (std::is_same_v<MATU, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>) {
      return std::make_pair(80, 4.690606234853613);
    } else
      static_assert(Dune::AlwaysFalse<MATU>::value,
                    "Expected results are not available for the given material and enhancedStrain type.");
  } else if constexpr (std::same_as<ES, EAS::DisplacementGradientTransposed>) {
    if constexpr (std::is_same_v<MATU, StVenantKirchhoff>) {
      return std::make_pair(80, 4.462882651578276);
    } else if constexpr (std::is_same_v<MATU, NeoHooke>) {
      return std::make_pair(80, 4.479501443200376);
    } else if constexpr (std::is_same_v<MATU, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>) {
      return std::make_pair(80, 4.683694983130406);
    } else
      static_assert(Dune::AlwaysFalse<MATU>::value,
                    "Expected results are not available for the given material and enhancedStrain type.");
  } else
    static_assert(Dune::AlwaysFalse<MATU>::value,
                  "Expected results are not available for the provided enhanced strain type.");
}

template <typename ES, typename MAT>
auto cantileverBeamResults3D(const MAT& mat) {
  using namespace Materials;
  if constexpr (std::same_as<ES, EAS::GreenLagrangeStrain>) {
    if constexpr (std::is_same_v<MAT, StVenantKirchhoff>) {
      return std::make_pair(80, 4.7692391315649365);
    } else if constexpr (std::is_same_v<MAT, NeoHooke>) {
      return std::make_pair(80, 4.781820664768682);
    } else if constexpr (std::is_same_v<MAT, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>) {
      return std::make_pair(80, 4.916861760967577);
    } else
      static_assert(Dune::AlwaysFalse<MAT>::value,
                    "Expected results are not available for the given material and enhancedStrain type.");
  } else if constexpr (std::same_as<ES, EAS::DisplacementGradient>) {
    if constexpr (std::is_same_v<MAT, StVenantKirchhoff>) {
      return std::make_pair(80, 4.751794459282368);
    } else if constexpr (std::is_same_v<MAT, NeoHooke>) {
      return std::make_pair(80, 4.763101490723167);
    } else if constexpr (std::is_same_v<MAT, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>) {
      return std::make_pair(80, 4.902880891417656);
    } else
      static_assert(Dune::AlwaysFalse<MAT>::value,
                    "Expected results are not available for the given material and enhancedStrain type.");
  } else if constexpr (std::same_as<ES, EAS::DisplacementGradientTransposed>) {
    if constexpr (std::is_same_v<MAT, StVenantKirchhoff>) {
      return std::make_pair(80, 4.745523824224497);
    } else if constexpr (std::is_same_v<MAT, NeoHooke>) {
      return std::make_pair(80, 4.755004891629581);
    } else if constexpr (std::is_same_v<MAT, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>) {
      return std::make_pair(80, 4.893238013810108);
    } else
      static_assert(Dune::AlwaysFalse<MAT>::value,
                    "Expected results are not available for the given material and enhancedStrain type.");
  } else
    static_assert(Dune::AlwaysFalse<MAT>::value,
                  "Expected results are not available for the provided enhanced strain type.");
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameter);
  Materials::NeoHooke matNH(matParameter);
  auto matBK       = Materials::makeBlatzKo(40.0);
  auto reducedMats = Dune::makeTupleVector(planeStrain(matSVK), planeStrain(matNH), planeStrain(matBK));
  auto materials   = Dune::makeTupleVector(matSVK, matNH, matBK);

  const auto testCantileverFunc = [&]<typename ES>() {
    Dune::Hybrid::forEach(reducedMats, [&t](const auto& mat) {
      t.subTest(cantileverBeamTest<2>(mat, skills(eas<ES>(4)), cantileverBeamResults<ES>(mat)));
    });

    auto params3d = std::same_as<ES, EAS::GreenLagrangeStrain> ? 21 : 9;
    Dune::Hybrid::forEach(materials, [&t, params3d](const auto& mat) {
      t.subTest(cantileverBeamTest<3>(mat, skills(eas<ES>(params3d)), cantileverBeamResults3D<ES>(mat)));
    });
  };

  testCantileverFunc.operator()<EAS::GreenLagrangeStrain>();
  testCantileverFunc.operator()<EAS::DisplacementGradient>();
  testCantileverFunc.operator()<EAS::DisplacementGradientTransposed>();

  return t.exit();
}
