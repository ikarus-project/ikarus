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
  } else if (elementType == "Q1P0") {
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
  } else {
    DUNE_THROW(Dune::NotImplemented, "Expected results are not implemented for the given element type");
  }
}

template <typename FEType>
auto testStressComputation(const FEType& fe, const auto& expectedStresses) {
  TestSuite t;

  auto n = fe.size();
  Eigen::VectorXd d;
  d.setZero(n);
  double lambda = 0.0;
  auto req      = typename FEType::Requirement(d, lambda);

  auto zeroStressState  = Eigen::Matrix<double, 2, 2>::Zero().eval();
  auto zeroStressStates = std::make_tuple(zeroStressState, zeroStressState, zeroStressState, zeroStressState);

  auto testStress = [&](const auto& d, const auto& expectedStresses) {
    const auto [expectedPK2Stress, expectedKirchhoffStress, expectedCauchyStress, expectedCauchyPolarStress] =
        expectedStresses;
    req.insertGlobalSolution(d);

    Dune::FieldVector<double, 2> pos{0.21132486540518711, 0.78867513459481287};
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
  const auto mat1 = Materials::makeOgden<1, Ikarus::PrincipalStretchTags::deviatoric>({matParameter.mu}, {2.0}, kappa,
                                                                                      Materials::VF1());
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
  t.subTest(testSingleElementImpl(basis1, skill2, gridView, "Q1E4"));
  t.subTest(testSingleElementImpl(basis1, skill3, gridView, "Q1H4"));
  t.subTest(testSingleElementImpl(basis1, skill4, gridView, "Q1HT4"));
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