// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"

#include <dune/common/parametertreeparser.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/fesettings.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/io/loadfesettings.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

auto createDummyElement() {
  auto grid     = createUGGridFromCorners<2>(CornerDistortionFlag::unDistorted, Dune::GeometryTypes::cube(2));
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis   = Ikarus::makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));
  auto element = gridView.begin<0>();

  using FEElement = Ikarus::LinearElastic<decltype(basis)>;
  return FEElement(basis, *element, 100, 0.2);
}

auto testFESettings() {
  TestSuite t("Test FE Settings");
  auto fe = createDummyElement();

  // Load PT
  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree("testfiles/fesettings.parset", parameterSet);

  auto settings  = Ikarus::Settings<decltype(fe)>(parameterSet.sub("FESettings"));
  auto container = settings.getContainer();

  t.check(container.nGP == 2) << "Setting is " << container.nGP;
  t.check(container.orderGP == 2) << "Setting is " << container.orderGP;

  constexpr auto key = refl::make_const_string("nGP");
  t.check(settings.get<key>() == container.nGP);

  fe.registerSettings(settings);

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testFESettings());

  return t.exit();
}
