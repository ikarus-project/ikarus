// SPDX-FileCopyrightText: 2021-2026 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "tests/src/testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/materialconversions.hh>
#include <ikarus/finiteelements/mechanics/materials/strainconversions.hh>
#include <ikarus/finiteelements/mechanics/materials/stressconversions.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;
using namespace Ikarus;

inline auto testMatrix() {
  return Eigen::Matrix3d{
      {0.600872, -0.179083,     0},
      {0.345678,  0.859121, -4e-5},
      { 0.00015,         0,     1}
  };
}

auto testMatPar() {
  double Emod = 1000;
  double nu   = 0.25;
  return YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
}

auto checkPK1AndPK2StressTensors() {
  using StrainTags::deformationGradient;

  TestSuite t("Check Two-Point Material Tensor: PK1 and PK2 tensor for SVK");
  spdlog::info("Testing " + t.name());

  auto matPar = testMatPar();
  auto svk    = Materials::StVenantKirchhoff(toLamesFirstParameterAndShearModulus(matPar));

  double epsilon                           = 1e-2;
  constexpr double scalingFactor           = 1e-2;
  constexpr double convergenceRateExpected = scalingFactor * scalingFactor;
  constexpr double tol                     = 1e-3;
  constexpr int numIncrements              = 5;

  std::vector<double> epsilons, errors;

  const Eigen::Matrix<double, 3, 3> F = testMatrix();
  Eigen::Matrix<double, 3, 3> deltaF;
  deltaF.setRandom();

  for (const auto counter : Dune::range(numIncrements)) {
    deltaF *= epsilon;
    const auto FPerturbed          = (F + deltaF).eval();
    const auto pk2Stress_exact     = svk.template stresses<deformationGradient, false>(F);
    const auto pk1Stress_exact     = transformStress<StressTags::PK2, StressTags::PK1>(pk2Stress_exact, F);
    const auto pk2Stress_perturbed = svk.template stresses<deformationGradient, false>(FPerturbed);
    const auto pk1Stress_perturbed = transformStress<StressTags::PK2, StressTags::PK1>(pk2Stress_perturbed, FPerturbed);
    const auto Cm                  = svk.template tangentModuli<deformationGradient, false>(F);
    const auto E                   = transformStrain<StrainTags::deformationGradient, StrainTags::greenLagrangian>(F);
    constexpr int dim              = 3;
    const auto A =
        transformTangentModuli<TangentModuliTags::Material, TangentModuliTags::TwoPoint>(Cm, pk2Stress_exact, F);

    Eigen::Matrix<double, 3, 3> pk1StressLin, pk2Stress;
    pk1StressLin.setZero();
    pk2Stress.setZero();

    for (const auto i : Dune::range(dim))
      for (const auto j : Dune::range(dim))
        for (const auto k : Dune::range(dim))
          for (const auto l : Dune::range(dim)) { // only works for SVK
            pk1StressLin(i, j) += A(i, j, k, l) * deltaF(k, l);
            pk2Stress(i, j) += Cm(i, j, k, l) * E(k, l);
          }

    Eigen::Matrix<double, 3, 3> pk1Error = pk1Stress_perturbed - pk1Stress_exact - pk1StressLin;
    double relErr                        = pk1Error.norm() / pk1Stress_perturbed.norm();

    std::cout << std::scientific << std::setprecision(8) << "epsilon " << epsilon << ", relErr " << relErr << std::endl;
    epsilons.push_back(epsilon);
    errors.push_back(relErr);
    checkApproxMatrices(t, pk2Stress_exact, pk2Stress, testLocation() + " Incorrect PK2 stress.", tol);
    epsilon = epsilon * 1e-2;
  }

  for (std::size_t i = 1; i < numIncrements; ++i) {
    t.check(Dune::FloatCmp::le(errors[i] / errors[i - 1], convergenceRateExpected + tol))
        << std::setprecision(16) << "Incorrect Scalar. Expected:\t" << convergenceRateExpected << " Actual:\t"
        << errors[i] / errors[i - 1] << ". The used tolerance was " << tol << " Incorrect error rate.";
  }

  return t;
}

auto checkTwoPointMaterialTensorSymmetry() {
  using StrainTags::deformationGradient;

  TestSuite t("Check Two-Point Material Tensor Implementation for Symmetry");
  spdlog::info("Testing " + t.name());

  auto matPar = testMatPar();
  auto nh     = Materials::NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

  constexpr double tol = 5e-14;

  auto F              = testMatrix();
  auto pk2Stress      = nh.template stresses<deformationGradient, false>(F);
  auto materialTensor = nh.template tangentModuli<deformationGradient, false>(F);
  constexpr int dim   = 3;
  const auto A =
      transformTangentModuli<TangentModuliTags::Material, TangentModuliTags::TwoPoint>(materialTensor, pk2Stress, F);

  for (const auto i : Dune::range(dim))
    for (const auto J : Dune::range(dim))
      for (const auto k : Dune::range(dim))
        for (const auto L : Dune::range(dim)) { // See Proposition 4.4 in Marsden and Hughes 1983 book
          const std::string& indexName = "i = " + std::to_string(i) + ", J = " + std::to_string(J) +
                                         ", k = " + std::to_string(k) + ", L = " + std::to_string(L);
          if (std::abs(A(i, J, k, L) - A(k, L, i, J)) > tol) { // due to floating-point issues
            checkScalars(t, A(i, J, k, L), A(k, L, i, J), " Incorrect entry for two-point tensor at " + indexName, tol);
          }
          checkScalars(t, materialTensor(i, J, k, L), materialTensor(J, i, k, L),
                       " Incorrect entry for material tensor at " + indexName, tol);
          checkScalars(t, materialTensor(i, J, k, L), materialTensor(i, J, L, k),
                       " Incorrect entry for material tensor at " + indexName, tol);
          checkScalars(t, materialTensor(i, J, k, L), materialTensor(k, L, i, J),
                       " Incorrect entry for material tensor at " + indexName, tol);
        }

  return t;
}

auto checkTwoPointMaterialTensorWithZeroStrains() {
  using StrainTags::deformationGradient;

  TestSuite t("Check Two-Point Material Tensor Implementation with Zero Strains");
  spdlog::info("Testing " + t.name());

  auto matPar = testMatPar();
  auto nh     = Materials::NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

  constexpr double tol = 1e-15;

  const Eigen::Matrix<double, 3, 3> F = Eigen::Matrix<double, 3, 3>::Identity();
  const auto pk2Stress                = nh.template stresses<deformationGradient, false>(F);
  const auto materialTensor           = nh.template tangentModuli<deformationGradient, false>(F);
  constexpr int dim                   = 3;
  const auto A =
      transformTangentModuli<TangentModuliTags::Material, TangentModuliTags::TwoPoint>(materialTensor, pk2Stress, F);

  Eigen::TensorFixedSize<double, Eigen::Sizes<dim, dim, dim, dim>> Aexp;
  Aexp.setZero();
  for (const auto i : Dune::range(dim))
    for (const auto J : Dune::range(dim))
      for (const auto k : Dune::range(dim))
        for (const auto L : Dune::range(dim)) {
          Aexp(i, J, k, L) += materialTensor(i, J, k, L); // true when F = I

          const std::string& indexName = "i = " + std::to_string(i) + ", J = " + std::to_string(J) +
                                         ", k = " + std::to_string(k) + ", L = " + std::to_string(L);
          checkScalars(t, A(i, J, k, L), Aexp(i, J, k, L), " Incorrect entry for two-point tensor at " + indexName,
                       tol);
        }

  return t;
}

auto checkTwoPointMaterialTensorWithDiagonalF(double lambda1, double lambda2, double lambda3) {
  using StrainTags::deformationGradient;

  TestSuite t("Check Two-Point Material Tensor Implementation for uniaxial tensile test (2D)");
  spdlog::info("Testing " + t.name() + " with lambda1 = " + std::to_string(lambda1) +
               " and lambda2 = " + std::to_string(lambda2) + " and lambda3 = " + std::to_string(lambda3));

  const auto matPar = toLamesFirstParameterAndShearModulus(testMatPar());
  const auto nh     = Materials::NeoHooke(matPar);
  const auto lambda = matPar.lambda;
  const auto mu     = matPar.mu;

  // A diagonal matrix chosen for F
  Eigen::Matrix<double, 3, 3> F;
  F.setZero();
  F(0, 0) = lambda1;
  F(1, 1) = lambda2;
  F(2, 2) = lambda3;

  constexpr double tol = 1e-15;
  constexpr int dim    = 3;

  const auto pk2Stress = nh.template stresses<deformationGradient, false>(F);
  const auto C         = nh.template tangentModuli<deformationGradient, false>(F);
  const auto A = transformTangentModuli<TangentModuliTags::Material, TangentModuliTags::TwoPoint>(C, pk2Stress, F);

  // A vector of index pairs vector the entries in tensors A or C are non-zero
  using Index4  = std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>;
  auto contains = [&](const std::vector<Index4>& v, const Index4& t) -> bool {
    return std::find(v.begin(), v.end(), t) != v.end();
  };
  std::vector<Index4> idxPairs;
  idxPairs.push_back({0, 0, 0, 0});
  idxPairs.push_back({0, 0, 1, 1});
  idxPairs.push_back({0, 0, 2, 2});
  idxPairs.push_back({0, 1, 0, 1});
  idxPairs.push_back({0, 1, 1, 0});
  idxPairs.push_back({0, 2, 0, 2});
  idxPairs.push_back({0, 2, 2, 0});
  idxPairs.push_back({1, 0, 0, 1});
  idxPairs.push_back({1, 0, 1, 0});
  idxPairs.push_back({1, 1, 0, 0});
  idxPairs.push_back({1, 1, 1, 1});
  idxPairs.push_back({1, 1, 2, 2});
  idxPairs.push_back({1, 2, 1, 2});
  idxPairs.push_back({1, 2, 2, 1});
  idxPairs.push_back({2, 0, 0, 2});
  idxPairs.push_back({2, 0, 2, 0});
  idxPairs.push_back({2, 1, 1, 2});
  idxPairs.push_back({2, 1, 2, 1});
  idxPairs.push_back({2, 2, 0, 0});
  idxPairs.push_back({2, 2, 1, 1});
  idxPairs.push_back({2, 2, 2, 2});

  const std::size_t numExpectedValues = idxPairs.size();

  std::vector<double> expectedAs, expectedCs;
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu + lambda + mu * lambda1 * lambda1) /
                       (lambda1 * lambda1));
  expectedAs.push_back((lambda) / (lambda1 * lambda2));
  expectedAs.push_back((lambda) / (lambda1 * lambda3));
  expectedAs.push_back(mu);
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda2));
  expectedAs.push_back(mu);
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda3));
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda2));
  expectedAs.push_back(mu);
  expectedAs.push_back((lambda) / (lambda1 * lambda2));
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu + lambda + mu * lambda2 * lambda2) /
                       (lambda2 * lambda2));
  expectedAs.push_back((lambda) / (lambda2 * lambda3));
  expectedAs.push_back(mu);
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda2 * lambda3));
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda3));
  expectedAs.push_back(mu);
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda2 * lambda3));
  expectedAs.push_back(mu);
  expectedAs.push_back((lambda) / (lambda1 * lambda3));
  expectedAs.push_back((lambda) / (lambda2 * lambda3));
  expectedAs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu + lambda + mu * lambda3 * lambda3) /
                       (lambda3 * lambda3));

  expectedCs.push_back((-2.0 * lambda * log(lambda1 * lambda2 * lambda3) + 2.0 * mu + lambda) /
                       (lambda1 * lambda1 * lambda1 * lambda1));
  expectedCs.push_back(lambda / (lambda1 * lambda1 * lambda2 * lambda2));
  expectedCs.push_back(lambda / (lambda1 * lambda1 * lambda3 * lambda3));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda1 * lambda2 * lambda2));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda1 * lambda2 * lambda2));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda1 * lambda3 * lambda3));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda1 * lambda3 * lambda3));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda1 * lambda2 * lambda2));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda1 * lambda2 * lambda2));
  expectedCs.push_back(lambda / (lambda1 * lambda1 * lambda2 * lambda2));
  expectedCs.push_back((-2.0 * lambda * log(lambda1 * lambda2 * lambda3) + 2.0 * mu + lambda) /
                       (lambda2 * lambda2 * lambda2 * lambda2));
  expectedCs.push_back(lambda / (lambda3 * lambda3 * lambda2 * lambda2));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda3 * lambda3 * lambda2 * lambda2));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda3 * lambda3 * lambda2 * lambda2));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda1 * lambda3 * lambda3));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda1 * lambda1 * lambda3 * lambda3));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda3 * lambda3 * lambda2 * lambda2));
  expectedCs.push_back((-lambda * log(lambda1 * lambda2 * lambda3) + mu) / (lambda3 * lambda3 * lambda2 * lambda2));
  expectedCs.push_back(lambda / (lambda1 * lambda1 * lambda3 * lambda3));
  expectedCs.push_back(lambda / (lambda3 * lambda3 * lambda2 * lambda2));
  expectedCs.push_back((-2.0 * lambda * log(lambda1 * lambda2 * lambda3) + 2.0 * mu + lambda) /
                       (lambda3 * lambda3 * lambda3 * lambda3));

  checkScalars(t, expectedAs.size(), numExpectedValues, " Incorrect size of expectedAs");
  checkScalars(t, expectedCs.size(), numExpectedValues, " Incorrect size of expectedCs");

  int counter = 0;
  for (const auto i : Dune::range(dim))
    for (const auto j : Dune::range(dim))
      for (const auto k : Dune::range(dim))
        for (const auto l : Dune::range(dim)) {
          Index4 currentPair           = {i, j, k, l};
          const std::string& indexName = "i = " + std::to_string(i) + ", j = " + std::to_string(j) +
                                         ", k = " + std::to_string(k) + ", l = " + std::to_string(l);
          if (contains(idxPairs, currentPair)) {
            checkScalars(t, A(i, j, k, l), expectedAs[counter], " Incorrect entry for two-point tensor at " + indexName,
                         tol);
            checkScalars(t, C(i, j, k, l), expectedCs[counter], " Incorrect entry for material tensor at " + indexName,
                         tol);
            counter++;
          } else {
            checkScalars(t, A(i, j, k, l), 0.0, " Incorrect entry for two-point tensor at " + indexName, tol);
            checkScalars(t, C(i, j, k, l), 0.0, " Incorrect entry for material tensor at " + indexName, tol);
          }
        }

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(checkPK1AndPK2StressTensors());
  t.subTest(checkTwoPointMaterialTensorSymmetry());
  t.subTest(checkTwoPointMaterialTensorWithZeroStrains());
  t.subTest(checkTwoPointMaterialTensorWithDiagonalF(1.15, 0.75, 1.39));
  t.subTest(checkTwoPointMaterialTensorWithDiagonalF(0.63, 1.95, 1.0));

  return t.exit();
}
