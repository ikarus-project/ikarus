// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/functionhelper.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <int order, int gridDim, typename ES>
auto easParametersFunc() {
  const std::string& errorMessage =
      "easParametersFunc is not implemented for the provided EnhancedStrainTag: " + ES::name();
  if constexpr (order == 1) {
    if constexpr (std::same_as<ES, EAS::GreenLagrangeStrain>) {
      std::array<int, (gridDim == 2) ? 4 : 3> easParameters;
      if constexpr (gridDim == 2)
        easParameters = {0, 4, 5, 7};
      else
        easParameters = {0, 9, 21};
      return easParameters;
    } else if constexpr (std::same_as<ES, EAS::DisplacementGradient> or
                         std::same_as<ES, EAS::DisplacementGradientTransposed>) {
      std::array<int, 2> easParameters;
      if constexpr (gridDim == 2)
        easParameters = {0, 4};
      else
        easParameters = {0, 9};
      return easParameters;
    } else
      DUNE_THROW(Dune::NotImplemented, errorMessage);
  } else {
    if constexpr (std::same_as<ES, EAS::GreenLagrangeStrain>) {
      static_assert(gridDim == 2);
      std::array<int, 2> easParameters = {0, 11};
      return easParameters;
    } else
      DUNE_THROW(Dune::NotImplemented, errorMessage);
  }
}

template <int order, int gridDim, typename ES, typename TestSuitType, typename MAT>
void easAutoDiffTest(TestSuitType& t, const MAT& mat) {
  if constexpr (order == 2 and not std::same_as<ES, EAS::GreenLagrangeStrain>)
    return;
  else {
    using namespace Dune::Functions::BasisFactory;
    const auto easParameters = easParametersFunc<order, gridDim, ES>();
    auto grid                = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
    auto gridView            = grid->leafGridView();

    for (const int numberOfInternalVariables : easParameters) {
      auto sk = skills(nonLinearElastic(mat), eas<ES>(numberOfInternalVariables));
      t.subTest(checkFESByAutoDiff(gridView, power<gridDim>(lagrange<order>()), sk,
                                   AffordanceCollections::elastoStatics,
                                   " (numberOfInternalVariables = " + std::to_string(numberOfInternalVariables) +
                                       ") and EnhancedStrainTag = " + ES::name()));
      if (numberOfInternalVariables != 0) {
        t.template checkThrow<Dune::NotImplemented>(
            [&]() {
              checkFESByAutoDiff(gridView, power<gridDim>(lagrange<3>()), sk, AffordanceCollections::elastoStatics,
                                 " (numberOfInternalVariables = " + std::to_string(numberOfInternalVariables) +
                                     ") and EnhancedStrainTag = " + ES::name());
            },
            "EAS method should have failed for a third-order basis.");
      }
    }
  }
}

template <int order, int gridDim, typename ES, typename MAT>
auto checkOrthogonalityCondition(const MAT& mat) {
  TestSuite t("Orthogonality Condition Test for EAS Method");
  if constexpr (order == 2 and not std::same_as<ES, EAS::GreenLagrangeStrain>)
    return t;
  else {
    using namespace Dune::Functions::BasisFactory;
    const auto easParameters  = easParametersFunc<order, gridDim, ES>();
    auto grid                 = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
    auto gridView             = grid->leafGridView();
    auto basis                = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<order>()));
    constexpr auto affordance = AffordanceCollections::elastoStatics;

    const auto orthogonalityFunc = [&]<typename FE>(const FE& fe) {
      const auto& localView  = fe.localView();
      const auto& element    = localView.element();
      auto rule              = Dune::QuadratureRules<double, gridDim>::rule(localView.element().type(), 2 * order);
      const auto& easVariant = fe.easVariant();
      const int numberOfInternalVariables = fe.numberOfInternalVariables();
      auto integrateFunctor               = [&]<typename EASF>(const EASF& easFunction) {
        typename EASF::AnsatzType integratedStrain;
        integratedStrain.setZero();
        for (const auto& gp : rule) {
          const auto M      = easFunction(gp.position());
          const double detJ = element.geometry().integrationElement(gp.position());
          if constexpr (std::same_as<ES, EAS::GreenLagrangeStrain>)
            integratedStrain += M * detJ * gp.weight();
          else
            for (const auto& H : M)
              integratedStrain += H * detJ * gp.weight();
        }
        t.check(integratedStrain.isZero()) << "Ansatz function for the enhanced strains in the EAS method should be "
                                                            "zero when integrated over the domain. (numberOfInternalVariables = " +
                                                  std::to_string(numberOfInternalVariables) +
                                                  ") and EnhancedStrainTag = " + ES::name();
      };
      easVariant(integrateFunctor);
    };

    for (const int numberOfInternalVariables : easParameters) {
      auto sk      = skills(nonLinearElastic(mat), eas<ES>(numberOfInternalVariables));
      auto fe      = makeFE(basis, sk);
      auto element = gridView.template begin<0>();
      fe.bind(*element);
      orthogonalityFunc(fe);
    }
    return t;
  }
}

template <int order, int gridDim, typename ES, typename MAT>
auto recoverNonlinearElastic(const MAT& mat) {
  TestSuite t(
      "Recover NonLinearElastic Test for EAS Method with numberOfInternalVariables = 0 with enhanced strain type as " +
      ES::name() + " and material type as " + Dune::className<MAT>());
  using namespace Dune::Functions::BasisFactory;
  auto grid                 = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
  auto gridView             = grid->leafGridView();
  auto basis                = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<order>()));
  constexpr auto affordance = AffordanceCollections::elastoStatics;
  auto element              = gridView.template begin<0>();
  auto nDOF                 = basis.flat().size();
  const double tol          = 1e-10;

  auto fe1 = makeFE(basis, skills(nonLinearElastic(mat)));
  auto fe2 = makeFE(basis, skills(nonLinearElastic(mat), eas<ES>(0)));
  fe1.bind(*element);
  fe2.bind(*element); // updateState need not be called as numberOfInternalVariables = 0

  Eigen::VectorXd d;
  d.setRandom(nDOF);
  d *= 0.01; // to avoid a negative determinant of deformation gradient
  double lambda = 0.0;
  auto req      = typename decltype(fe1)::Requirement(d, lambda);

  Eigen::MatrixXd K1, K2;
  K1.setZero(nDOF, nDOF);
  K2.setZero(nDOF, nDOF);

  Eigen::VectorXd R1, R2;
  R1.setZero(nDOF);
  R2.setZero(nDOF);

  calculateMatrix(fe1, req, affordance.matrixAffordance(), K1);
  calculateMatrix(fe2, req, affordance.matrixAffordance(), K2);

  calculateVector(fe1, req, affordance.vectorAffordance(), R1);
  calculateVector(fe2, req, affordance.vectorAffordance(), R2);

  checkApproxMatrices(t, K1, K2, testLocation() + "\nIncorrect stiffness matrices.", tol);
  checkSymmetricMatrix(t, K1, tol, "K1");
  checkSymmetricMatrix(t, K2, tol, "K2");
  checkApproxVectors(t, R1, R2, testLocation() + "\nIncorrect residual vectors.", tol);
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t("Nonlinear EAS Test");
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameter);
  Materials::NeoHooke matNH(matParameter);
  auto matBK         = Materials::makeBlatzKo(40.0);
  auto reducedMatSVK = planeStrain(matSVK);
  auto reducedMatNH  = planeStrain(matNH);
  auto reducedMatBK  = planeStrain(matBK);

  const auto testNonlinearEASFunc = [&]<typename ES>() {
    easAutoDiffTest<1, 2, ES>(t, reducedMatSVK);
    easAutoDiffTest<1, 3, ES>(t, matSVK);
    easAutoDiffTest<1, 2, ES>(t, reducedMatNH);
    easAutoDiffTest<1, 3, ES>(t, matNH);
    easAutoDiffTest<1, 2, ES>(t, reducedMatBK);
    easAutoDiffTest<1, 3, ES>(t, matBK);
    easAutoDiffTest<2, 2, ES>(t, reducedMatSVK);
    easAutoDiffTest<2, 2, ES>(t, reducedMatNH);
    easAutoDiffTest<2, 2, ES>(t, reducedMatBK);

    t.subTest(recoverNonlinearElastic<1, 2, ES>(reducedMatSVK));
    t.subTest(recoverNonlinearElastic<1, 3, ES>(matSVK));
    t.subTest(recoverNonlinearElastic<1, 2, ES>(reducedMatNH));
    t.subTest(recoverNonlinearElastic<1, 3, ES>(matNH));
    t.subTest(recoverNonlinearElastic<2, 2, ES>(reducedMatSVK));
    t.subTest(recoverNonlinearElastic<2, 2, ES>(reducedMatNH));

    t.subTest(checkOrthogonalityCondition<1, 2, ES>(reducedMatSVK));
    t.subTest(checkOrthogonalityCondition<1, 3, ES>(matSVK));
    t.subTest(checkOrthogonalityCondition<2, 2, ES>(reducedMatSVK));
  };

  testNonlinearEASFunc.operator()<EAS::GreenLagrangeStrain>();
  testNonlinearEASFunc.operator()<EAS::DisplacementGradient>();
  testNonlinearEASFunc.operator()<EAS::DisplacementGradientTransposed>();

  return t.exit();
}