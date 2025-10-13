// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/displacementpressure.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolverfactory.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using Dune::TestSuite;

namespace Testing {
template <typename FEType>
auto testEigenValuesQ1P0(const FEType& fe) {
  TestSuite t;

  auto n = fe.size();

  Eigen::VectorXd d;
  d.setZero(n);
  double lambda = 0.0;
  auto req      = typename FEType::Requirement(d, lambda);

  auto testEVs = [&](const auto& d, const auto& evsExpected, const auto& normfexpected, const auto& expectedStress) {
    req.insertGlobalSolution(d);
    Eigen::MatrixXd k;
    k.setZero(n, n);

    Eigen::VectorXd f;
    f.setZero(n);

    calculateMatrix(fe, req, Ikarus::MatrixAffordance::stiffness, k);
    calculateVector(fe, req, Ikarus::VectorAffordance::forces, f);

    auto essaK = makeIdentitySymEigenSolver<EigenValueSolverType::Eigen>(k);
    essaK.compute();
    auto eigenValuesComputed = essaK.eigenvalues();

    checkApproxVectors(t, eigenValuesComputed, evsExpected, " Incorrect eigenvalues of K ", 1e-10);
    checkScalars(t, f.norm(), normfexpected, " Incorrect norm of f", 1e-10);

    Dune::FieldVector<double, 2> firstQP{0.2113248654052, 0.21132486540520};
    auto S = fe.template calculateAt<ResultTypes::PK2Stress>(req, firstQP);

    checkApproxMatrices(t, S.asMat(), expectedStress, "Incorrect Stress", 1e-8);
  };

  auto expectedEigenValues0 = Eigen::Vector<double, 9>{
      -0.00899964037263862, -2.1729667857278618e-14, 2.946341739174961e-14, 1.1436002929121896e-13, 223.72258617281094,
      260.9992542878448,    260.9992542878448,       671.1409395973144,     671.1409395973158};
  auto expectedEigenValuesD = Eigen::Vector<double, 9>{
      -35.283103786777886, -0.01060548937017789, -5.267913592314902e-15, 4.55693342945812e-14, 142.21651420720613,
      198.0477505746059,   382.0279746326626,    715.9234881273067,      1102.7590489808263};

  auto expectedStress0 = Eigen::Matrix<double, 2, 2>::Zero().eval();
  auto expectedStressD = Eigen::Matrix<double, 2, 2>{
      {79.83415515,   43.14730543},
      {43.14730543, -115.61657004}
  };

  testEVs(d, expectedEigenValues0, 8.038873388460925e-14, expectedStress0);

  d << 0.0, 0.1, 0.2, 0.3, 0.0, 0.1, 0.4, -0.1, 1e-5;
  testEVs(d, expectedEigenValuesD, 278.35552561152656, expectedStressD);

  return t;
};

template <typename FEType>
auto testEigenValuesQ2P1(const FEType& fe) {
  TestSuite t;

  auto n = fe.size();

  Eigen::VectorXd d;
  d.setZero(n);
  double lambda = 0.0;
  auto req      = typename FEType::Requirement(d, lambda);

  auto testEVs = [&](const auto& d, const auto& evsExpected, const auto& normfexpected) {
    req.insertGlobalSolution(d);
    Eigen::MatrixXd k;
    k.setZero(n, n);

    Eigen::VectorXd f;
    f.setZero(n);

    calculateMatrix(fe, req, Ikarus::MatrixAffordance::stiffness, k);
    calculateVector(fe, req, Ikarus::VectorAffordance::forces, f);

    auto essaK = makeIdentitySymEigenSolver<EigenValueSolverType::Eigen>(k);
    essaK.compute();
    auto eigenValuesComputed = essaK.eigenvalues();

    checkApproxVectors(t, eigenValuesComputed, evsExpected, " Incorrect eigenvalues of K ", 1e-10);
    checkScalars(t, f.norm(), normfexpected, " Incorrect norm of f", 1e-10);
  };

  auto expectedEigenValues0 = Eigen::Vector<double, 22>{
      -0.0022499662832745166,  -0.0007499953170165504, -0.0007499953170165504, -0.00019269216525575703,
      -1.2675836675755374e-13, 1.7663987885579953e-14, 1.7663987885579953e-14, 74.57121551081288,
      95.04370854203668,       95.04370854203673,      137.4532366264664,      202.61599074588514,
      340.8307663627087,       340.83076636270874,     715.8836689038023,      856.2954605331539,
      984.3400447427299,       1100.4311758193091,     1177.6687024698504,     1177.6687024698533,
      2614.645487083809,       2614.6454870838124};
  auto expectedEigenValuesD = Eigen::Vector<double, 22>{
      -4427.013285300488,     -0.002040646174187323,  -0.0008330950847293147, -0.00039650220441031107,
      -3.8373673614607304e-5, 1.0528132510634611e-12, 3.646346256005411e-12,  10.396215358782568,
      78.67635569586426,      112.67502157113005,     168.1386004850153,      283.156382121772,
      334.1004165735087,      673.497820795596,       693.4161372086264,      1013.8903723881498,
      1089.7042768037954,     1379.5488153262536,     2097.470724366695,      2354.0018864800427,
      7806.442489077462,      85336.8479396152};

  testEVs(d, expectedEigenValues0, 8.038873388460929e-14);

  return t;
};

template <bool continuous = true>
auto eigenvaluesForAssemblerTest() {
  if constexpr (not continuous) {
    auto expectedEigenValues0 = Eigen::Vector<double, 66>{
        -0.0005624956097084438,  -0.0004951385262628075,  -0.0004951385262461518,  -0.00042785067620338313,
        -0.0001784246342901316,  -0.00014783016412167144, -0.00014783016412167144, -0.00010974263446594822,
        -0.00010653812258366693, -7.205523162150919e-5,   -7.205523152250576e-5,   -5.796738527550611e-5,
        -2.3878371505434913e-5,  -2.387837142053355e-5,   -2.200639820133884e-5,   -1.7814161293814225e-5,
        8.182344433071662e-14,   1.231529015048884e-13,   1.231529015048884e-13,   41.78277474873721,
        41.78277474873722,       42.425032642478676,      68.58772307883274,       103.18922154911795,
        106.32542331746008,      126.27886231078458,      126.27886231078554,      207.29073610476678,
        221.51812600114903,      221.51812600114917,      223.46972629254708,      283.21195067695896,
        344.0908281044,          427.79494208839,         484.7446122592129,       484.74461225921374,
        619.4728850760116,       619.4728850760123,       697.5631959027543,       697.5631959027555,
        747.8764214411284,       753.7244704957925,       901.5345916351583,       992.6592947122655,
        1057.5000851318523,      1057.5000851318553,      1121.4024147958544,      1131.1845498575553,
        1131.1845498575622,      1166.687772067528,       1266.6885048475592,      1413.9026718150621,
        1437.1500213593963,      1437.150021359399,       1489.6128625178426,      1594.466522102227,
        1856.8509427099539,      1856.8509427099577,      2455.1718843709,         2455.1718843709023,
        2582.4002815506465,      2582.8808631531842,      2818.7318410403245,      2842.578569101952,
        2885.9456063278435,      2885.945606327856};

    auto expectedEigenValuesD = Eigen::Vector<double, 66>{
        -37239.96585874105,      -0.0005581058118415022,  -0.0004952459762606795,  -0.0004628869304182488,
        -0.0003302214097235077,  -0.00017667090408643368, -0.00015286961708233158, -0.00013660386182210818,
        -0.00011454998561127305, -8.963363060498917e-5,   -7.294189114918169e-5,   -6.164162821388082e-5,
        -3.188034220709957e-5,   -2.432748707400558e-5,   -2.2437776814135436e-5,  -1.93064689135687e-5,
        -4.520521781258805e-6,   -3.720210793442427e-11,  1.748745366599946e-12,   0.7516464424101633,
        42.176270956177845,      42.98024989354218,       69.28379195655454,       71.18231645917474,
        108.97947229323614,      118.73119728427093,      159.88970027609446,      204.0024568958243,
        215.9157227059967,       229.15767575095853,      290.4493812713574,       310.9622951465548,
        353.32692859434064,      467.26914346210685,      503.905448263122,        576.7274374293804,
        595.3260062242558,       683.5361836921437,       694.9442492803614,       724.124244737345,
        755.2165676276632,       910.1405765704194,       965.6288628894806,       1032.4627260793807,
        1059.3359129809623,      1100.6119186614137,      1139.3678124380538,      1159.2923019431346,
        1258.9715879036842,      1377.424802621225,       1403.4370963722913,      1435.8601530832302,
        1490.6980897887152,      1576.7283810458478,      1836.466165393944,       1859.339455640058,
        2463.3557921492184,      2467.651389099999,       2514.5143682732783,      2598.843107043307,
        2739.829815030995,       2797.5516579402033,      2865.6195269656214,      2875.9970125505397,
        42421.33981829227,       2.727112499985963e6};

    return std::make_pair(expectedEigenValues0, expectedEigenValuesD);
  } else {
    auto expectedEigenValues0 = Eigen::Vector<double, 59>{
        -0.0013117301430375363, -0.0005860444513958874, -0.000586044451308379,   -0.00035684818922466485,
        -0.0002458414083269222, -0.0001480623747877602, -0.00013683869333266805, -0.00013683869321414455,
        -5.432120523508456e-5,  -1.791765953808241e-13, -5.0128880983085575e-14, 8.351667685129325e-14,
        41.7829114859163,       41.78291148591643,      42.425157675845064,      68.58841522503216,
        103.18913903196754,     106.32535360291087,     126.27876966708668,      126.2787696670867,
        207.2907923212113,      221.51812069199627,     221.51812069199653,      223.46973517536483,
        283.21187441973217,     344.0908335947386,      427.79493372659675,      484.7446090770773,
        484.7446090770778,      619.4728842158041,      619.4728842158047,       697.5631924600625,
        697.563192460065,       747.876427691951,       753.7244657435759,       901.5345891484162,
        992.6593062786966,      1057.500085684075,      1057.500085684077,       1121.4024170431658,
        1131.1845493021913,     1131.1845493021935,     1166.6877712676362,      1266.688503048387,
        1413.9026911556869,     1437.1500202533,        1437.1500202533027,      1489.6128035887828,
        1594.4665615663098,     1856.8509265354414,     1856.850926535448,       2455.1718741137593,
        2455.1718741137624,     2582.400297729741,      2582.880873499452,       2818.731811300359,
        2842.578541926258,      2885.9455882169595,     2885.945588216972};

    auto expectedEigenValuesD = Eigen::Vector<double, 59>{
        -37240.025763563404,     -0.0012951996293419316,  -0.00059992332610356,    -0.0005318951210979956,
        -0.00035723855584915003, -0.00021423360150700044, -0.00014473979440250608, -0.00013725700070047232,
        -6.482132084577299e-5,   -2.191802822616197e-5,   -4.999614849556142e-11,  3.644851765145945e-12,
        0.7440821392179334,      42.17246948329782,       42.96014201936403,       69.25889939153329,
        71.17117762743096,       108.97330302551889,      118.72312942558649,      159.86930936215464,
        203.9818103964321,       215.90829054060652,      229.1487468533353,       290.44758783004767,
        310.9339671429349,       353.31367260181133,      467.26922388421514,      503.889136409957,
        576.7404778206267,       595.3296465889384,       683.5401531138409,       694.9402004891882,
        724.1366894559294,       755.2419270686445,       910.1507116108429,       965.6388661830865,
        1032.4734104441961,      1059.3574234063753,      1100.618085100647,       1139.3999997039477,
        1159.2952715765039,      1258.989226287927,       1377.4142188914857,      1403.4587590355056,
        1435.8648423012532,      1490.692575486046,       1576.7393895097935,      1836.4602171363522,
        1859.3454957676674,      2463.360476712742,       2467.6519686231463,      2514.5131382059308,
        2598.853293908254,       2739.8304555252157,      2797.5432595213583,      2865.6182365648187,
        2876.0025388311165,      42421.40012247569,       2.72711250119609e6};

    return std::make_pair(expectedEigenValues0, expectedEigenValuesD);
  }
}
} // namespace Testing

template <int pD, int pP, bool continuous = true>
auto testSingleElement() {
  TestSuite t("Eigenvalue test: u-p Element");
  spdlog::info("Testing " + t.name() + " with pD = " + std::to_string(pD) + " and pP = " + std::to_string(pP));

  // SETUP
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 1000.0, .nu = 0.49});
  auto kappa        = convertLameConstants(matParameter).toBulkModulus();
  auto mat          = Materials::makeOgden<1, Ikarus::PrincipalStretchTags::deviatoric>({matParameter.mu}, {2.0}, kappa,
                                                                                        Materials::VF12());
  auto reducedMat   = planeStrain(mat);
  auto grid         = createGrid<Grids::Yasp>(1, 1);
  auto gridView     = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = [&]() {
    if constexpr (continuous)
      return Ikarus::makeBasis(
          gridView, composite(power<2>(lagrange<pD>(), FlatInterleaved{}), lagrange<pP>(), BlockedLexicographic{}));
    else
      return Ikarus::makeBasis(
          gridView,
          composite(power<2>(lagrange<pD>(), FlatInterleaved{}), lagrangeDG<pP>(),
                    BlockedLexicographic{})); // for a single element, lagrange and lagrangeDG should be the same
  }();
  auto sk      = skills(displacementPressure(reducedMat));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  auto& fe = fes.front();

  auto nDOF = basis.flat().size();
  auto n    = fe.size();

  t.check(n == nDOF);

  Eigen::VectorXd d;
  d.setZero(n);

  // TESTS
  if (pD == 1) {
    t.subTest(Testing::testEigenValuesQ1P0(fe));
  } else {
    t.subTest(Testing::testEigenValuesQ2P1(fe));
  }
  return t;
}

template <int pD, int pP, bool continuous = true>
auto testAssembler() {
  TestSuite t("Assembler Test");
  spdlog::info("Testing " + t.name() + " with pD = " + std::to_string(pD) + " and pP = " + std::to_string(pP) +
               " with continuous being " + std::to_string(continuous));

  // SETUP
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 1000.0, .nu = 0.49});
  auto kappa        = convertLameConstants(matParameter).toBulkModulus();
  auto mat          = Materials::makeOgden<1, Ikarus::PrincipalStretchTags::deviatoric>({matParameter.mu}, {2.0}, kappa,
                                                                                        Materials::VF12());
  auto reducedMat   = planeStrain(mat);
  auto grid         = createGrid<Grids::Yasp>(2, 2);
  auto gridView     = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = [&]() {
    if constexpr (continuous)
      return Ikarus::makeBasis(
          gridView, composite(power<2>(lagrange<pD>(), FlatInterleaved{}), lagrange<pP>(), BlockedLexicographic{}));
    else
      return Ikarus::makeBasis(
          gridView, composite(power<2>(lagrange<pD>(), FlatInterleaved{}), lagrangeDG<pP>(), BlockedLexicographic{}));
  }();
  auto sk      = skills(displacementPressure(reducedMat));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  Ikarus::DirichletValues dirichletValues(basis.flat());
  auto sparseFlatAssembler = makeSparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;
  auto req      = typename FEType::Requirement(d, lambda);

  sparseFlatAssembler->bind(req, Ikarus::AffordanceCollections::elastoStatics, DBCOption::Full);

  auto K = sparseFlatAssembler->matrix().toDense();

  auto essaK = makeIdentitySymEigenSolver<EigenValueSolverType::Eigen>(K);
  essaK.compute();
  auto eigenValuesComputed = essaK.eigenvalues();

  const auto [expectedEigenValues0, expectedEigenValuesD] = Testing::eigenvaluesForAssemblerTest<continuous>();

  checkApproxVectors(t, eigenValuesComputed, expectedEigenValues0, " Incorrect eigenvalues of K ", 1e-10);

  // "Slightly deformed" state
  d[0] = 0.3; // First d dof (first element)
  if constexpr (not continuous) {
    d[50] = 0.1; // first p dof (first element)
    d[54] = 0.3; // First p dof (second element)
  }

  auto reqD = typename FEType::Requirement(d, lambda);
  sparseFlatAssembler->bind(reqD, Ikarus::AffordanceCollections::elastoStatics, DBCOption::Full);

  K = sparseFlatAssembler->matrix().toDense();

  auto essaKD = makeIdentitySymEigenSolver<EigenValueSolverType::Eigen>(K);
  essaKD.compute();
  auto eigenValuesComputedD = essaKD.eigenvalues();

  checkApproxVectors(t, eigenValuesComputedD, expectedEigenValuesD, " Incorrect eigenvalues of K ", 1e-10);

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;

  t.subTest(testSingleElement<1, 0>());
  t.subTest(testSingleElement<2, 1>());
  t.subTest(testSingleElement<2, 1, false>());
  t.subTest(testAssembler<2, 1>());
  t.subTest(testAssembler<2, 1, false>());

  return t.exit();
}