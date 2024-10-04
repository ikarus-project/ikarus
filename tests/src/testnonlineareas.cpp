// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "testcommon.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <int gridDim, typename TestSuitType, typename MAT>
void easAutoDiffTest(const TestSuitType& t, const MAT& mat) {
  using namespace Dune::Functions::BasisFactory;
  std::array<int, (gridDim == 2) ? 4 : 3> easParameters;
  if constexpr (gridDim == 2)
    easParameters = {0, 4, 5, 7};
  else
    easParameters = {0, 9, 21};

  auto grid = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
  grid->globalRefine(2);
  auto gridView = grid->leafGridView();

  for (const int numberOfEASParameters : easParameters) {
    auto sk = skills(nonLinearElastic(mat), eas(numberOfEASParameters));
    checkFESByAutoDiff(gridView, power<gridDim>(lagrange<1>()), sk, AffordanceCollections::elastoStatics,
                       " (numberOfEASParameters = " + std::to_string(numberOfEASParameters) + ")");
  }
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t("Nonlinear EAS Test");
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  StVenantKirchhoff matSVK(matParameter);
  auto reducedMat = planeStress(matSVK);
  easAutoDiffTest<2>(t, reducedMat);
  easAutoDiffTest<3>(t, matSVK);
  return t.exit();
}