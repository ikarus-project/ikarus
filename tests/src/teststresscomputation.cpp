// SPDX-FileCopyrightText: 2021-2026 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/assumedstress.hh>
#include <ikarus/finiteelements/mechanics/displacementpressure.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/io/resultevaluators.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using Dune::TestSuite;

namespace Testing {

auto expectedResults(const std::string& elementType) {
  using MatrixType = Eigen::Matrix<double, 2, 2>;
  if (elementType == "Q1") {
    auto expectedPK2Stress = MatrixType{
        {3445.9817205200880, 53.091836258543370},
        {53.091836258543370, 3205.4837845288780}
    };
    auto expectedKirchhoffStress = MatrixType{
        {5329.2168543939660, 678.97252282322170},
        {678.97252282322170, 2743.6402549203340}
    };
    auto expectedCauchyStress = MatrixType{
        {4706.2253630377540, 599.59986523756390},
        {599.59986523756390, 2422.9055989927070}
    };
    auto expectedCauchyPolarStress = MatrixType{
        { 4164.1653462527930, -1141.6598820225230},
        {-1141.6598820225230,  2964.9656157776670}
    };
    return std::make_tuple(expectedPK2Stress, expectedKirchhoffStress, expectedCauchyStress, expectedCauchyPolarStress);
  } else if (elementType == "Q1E4") {
    auto expectedPK2Stress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedKirchhoffStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedCauchyStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedCauchyPolarStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    return std::make_tuple(expectedPK2Stress, expectedKirchhoffStress, expectedCauchyStress, expectedCauchyPolarStress);
  } else if (elementType == "Q1H4") {
    auto expectedPK2Stress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedKirchhoffStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedCauchyStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedCauchyPolarStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    return std::make_tuple(expectedPK2Stress, expectedKirchhoffStress, expectedCauchyStress, expectedCauchyPolarStress);
  } else if (elementType == "Q1HT4") {
    auto expectedPK2Stress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedKirchhoffStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedCauchyStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    auto expectedCauchyPolarStress = MatrixType{
        {0.0, 0.0},
        {0.0, 0.0}
    };
    return std::make_tuple(expectedPK2Stress, expectedKirchhoffStress, expectedCauchyStress, expectedCauchyPolarStress);
  } else if (elementType == "Q1S5") {
    auto expectedPK2Stress = MatrixType{
        {1769.6604462179770, 33.557046979865760},
        {33.557046979865760, 1662.0213894432550}
    };
    auto expectedKirchhoffStress = MatrixType{
        {2737.4720926549360, 356.48208114234320},
        {356.48208114234320, 1423.6059549829540}
    };
    auto expectedCauchyStress = MatrixType{
        {2417.4585018882210, 314.80892175692160},
        {314.80892175692160, 1257.1848050785050}
    };
    auto expectedCauchyPolarStress = MatrixType{
        {2152.1305752402840000, -580.1368484048580},
        {   -580.1368484048580, 1522.5127317264410}
    };
    return std::make_tuple(expectedPK2Stress, expectedKirchhoffStress, expectedCauchyStress, expectedCauchyPolarStress);
  } else if (elementType == "Q1P0") {
    auto expectedPK2Stress = MatrixType{
        {79.83415515,   43.14730543},
        {43.14730543, -115.61657004}
    };
    auto expectedKirchhoffStress = MatrixType{
        {127.52615425356670,  56.258411673211170},
        {56.258411673211170, -86.709956697111760}
    };
    auto expectedCauchyStress = MatrixType{
        {112.61820226061170,  49.681739575372270},
        {49.681739575372270, -76.573464466808380}
    };
    auto expectedCauchyPolarStress = MatrixType{
        { 67.704108472273930, -94.595833363710030},
        {-94.595833363710030,  -31.65937067847060}
    };
    return std::make_tuple(expectedPK2Stress, expectedKirchhoffStress, expectedCauchyStress, expectedCauchyPolarStress);
  } else {
    DUNE_THROW(Dune::NotImplemented, "Expected results are not implemented for the given element type");
  }
}

template <typename FEType>
auto testStressComputation(FEType& fe, const auto& expectedStresses) {
  TestSuite t;

  auto n = fe.size();
  Eigen::VectorXd d;
  d.setZero(n);
  double lambda = 0.0;
  auto req      = typename FEType::Requirement(d, lambda);

  auto zeroStressState  = Eigen::Matrix<double, 2, 2>::Zero().eval();
  auto zeroStressStates = std::make_tuple(zeroStressState, zeroStressState, zeroStressState, zeroStressState);

  auto testStress = [&](const auto& disp, const auto& expectedStresses) {
    const auto [expectedPK2Stress, expectedKirchhoffStress, expectedCauchyStress, expectedCauchyPolarStress] =
        expectedStresses;

    // In a nonlinear analysis, the solver takes care of updating the state when CORRECTION_UPDATED is broadcasted. In
    // other words, Rtilde at current iteration with old displacements is used to update internal variables (alpha/beta)
    // and only after that the displacements are updated. But here, it has to be done manually and this goal is
    // acheieved by manually updating the requirements later.

    fe.updateState(req, disp); // here disp = correction vector (DeltaD) as original displacements are zero
    req.insertGlobalSolution(disp);

    Dune::FieldVector<double, 2> pos{0.21132486540518711, 0.21132486540518711};
    const auto S                              = fe.template calculateAt<ResultTypes::PK2Stress>(req, pos);
    const auto tau                            = fe.template calculateAt<ResultTypes::kirchhoffStress>(req, pos);
    const auto sigma                          = fe.template calculateAt<ResultTypes::cauchyStress>(req, pos);
    const Eigen::Vector<double, 3> sigmaVoigt = sigma.asVec();

    auto polarStress = ResultEvaluators::PolarStress({0.0, 0.0});
    Eigen::Matrix<double, 2, 2> sigmaPolar;
    sigmaPolar.setZero();
    sigmaPolar << polarStress(sigmaVoigt, pos, fe, 0), polarStress(sigmaVoigt, pos, fe, 2),
        polarStress(sigmaVoigt, pos, fe, 2), polarStress(sigmaVoigt, pos, fe, 1);

    const double tol = 1e-8;

    checkApproxMatrices(t, S.asMat(), expectedPK2Stress, "Incorrect PK2Stress", tol);
    checkApproxMatrices(t, tau.asMat(), expectedKirchhoffStress, "Incorrect kirchhoffStress", tol);
    checkApproxMatrices(t, sigma.asMat(), expectedCauchyStress, "Incorrect cauchyStress", tol);
    checkApproxMatrices(t, sigmaPolar, expectedCauchyPolarStress, "Incorrect polarCauchyStress", tol);

    checkSymmetricMatrix(t, S.asMat(), tol, "PK2Stress not symmetric");
    checkSymmetricMatrix(t, tau.asMat(), tol, "kirchhoffStress not symmetric");
    checkSymmetricMatrix(t, sigma.asMat(), tol, "cauchyStress not symmetric");
    checkSymmetricMatrix(t, sigmaPolar, tol, "polarCauchyStress not symmetric");
  };

  testStress(d, zeroStressStates);

  if (n == 9)
    d << 0.0, 0.1, 0.2, 0.3, 0.0, 0.1, 0.4, -0.1, 1e-5;
  else
    d << 0.0, 0.1, 0.2, 0.3, 0.0, 0.1, 0.4, -0.1;

  testStress(d, expectedStresses);

  return t;
}
} // namespace Testing

template <typename B, typename S, typename GV>
auto testSingleElementImpl(const B& basis, const S& skill, const GV& gridView, const std::string& elementType = "") {
  TestSuite t("Stress Computation Test with element type " + elementType);
  spdlog::info("Testing " + t.name());

  using FEType = decltype(makeFE(basis, skill));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, skill));
    fes.back().bind(ge);
  }

  auto& fe = fes.front();

  auto nDOF = basis.flat().size();
  auto n    = fe.size();

  t.check(n == nDOF);
  t.subTest(Testing::testStressComputation(fe, Testing::expectedResults(elementType)));

  return t;
}

auto testSingleElement() {
  TestSuite t;

  auto grid     = createGrid<Grids::Yasp>(1, 1);
  auto gridView = grid->leafGridView();

  const auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 1000.0, .nu = 0.49});
  const auto kappa        = convertLameConstants(matParameter).toBulkModulus();
  const auto mat1         = Materials::StVenantKirchhoff(matParameter);
  const auto mat2 = Materials::makeOgden<1, Ikarus::PrincipalStretchTags::deviatoric>({matParameter.mu}, {2.0}, kappa,
                                                                                      Materials::VF12());
  const auto reducedMat1 = planeStrain(mat1);
  const auto reducedMat2 = planeStrain(mat2);

  using namespace Dune::Functions::BasisFactory;
  const auto basis1 = Ikarus::makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved{}));
  const auto basis2 = Ikarus::makeBasis(
      gridView, composite(power<2>(lagrange<1>(), FlatInterleaved{}), lagrange<0>(), BlockedLexicographic{}));

  const auto skill1 = skills(nonLinearElastic(reducedMat1));
  const auto skill2 = skills(nonLinearElastic(reducedMat1), eas<EAS::GreenLagrangeStrain>(4));
  const auto skill3 = skills(nonLinearElastic(reducedMat1), eas<EAS::DisplacementGradient>(4));
  const auto skill4 = skills(nonLinearElastic(reducedMat1), eas<EAS::DisplacementGradientTransposed>(4));
  const auto skill5 = skills(nonLinearElastic(reducedMat1), assumedStress<PS::PK2Stress>(5));
  const auto skill6 = skills(displacementPressure(reducedMat2));

  t.subTest(testSingleElementImpl(basis1, skill1, gridView, "Q1"));
  // t.subTest(testSingleElementImpl(basis1, skill2, gridView, "Q1E4"));
  // t.subTest(testSingleElementImpl(basis1, skill3, gridView, "Q1H4"));
  // t.subTest(testSingleElementImpl(basis1, skill4, gridView, "Q1HT4"));
  t.subTest(testSingleElementImpl(basis1, skill5, gridView, "Q1S5"));
  t.subTest(testSingleElementImpl(basis2, skill6, gridView, "Q1P0"));
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;

  t.subTest(testSingleElement());

  return t.exit();
}