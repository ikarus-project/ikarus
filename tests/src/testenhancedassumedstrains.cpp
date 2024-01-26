// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "resultcollection.hh"
#include "testeas.hh"
#include "testfeelement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template <typename Basis>
using EASElement = Ikarus::EnhancedAssumedStrains<Ikarus::LinearElastic<Basis>>;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("EAS");

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis = power<2>(lagrange<1>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis = power<3>(lagrange<1>(), FlatInterleaved());
  constexpr auto randomlyDistorted      = CornerDistortionFlag::randomlyDistorted;
  constexpr auto unDistorted            = CornerDistortionFlag::unDistorted;

  t.subTest(testFEElement<EASElement>(firstOrderLagrangePrePower2Basis, "EAS", randomlyDistorted,
                                      Dune::ReferenceElements<double, 2>::cube(), checkJacobianFunctor));

  t.subTest(testFEElement<EASElement>(firstOrderLagrangePrePower3Basis, "EAS", randomlyDistorted,
                                      Dune::ReferenceElements<double, 3>::cube(), checkJacobianFunctor));
  t.subTest(testFEElement<EASElement>(
      firstOrderLagrangePrePower2Basis, "EAS", unDistorted, Dune::ReferenceElements<double, 2>::cube(),
      checkCalculateAtFunctorFactory<Ikarus::ResultType::linearStress>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<Ikarus::ResultType::linearStress>()));

  return t.exit();
}
