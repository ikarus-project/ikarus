// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "tests/src/testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/materialconversions.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;
using namespace Ikarus;

inline auto testMatrix1() {
  return Eigen::Matrix3d{
      {0.600872, -0.179083,     0},
      {0.345678,  0.859121, -4e-5},
      { 0.00015,         0,     1}
  };
}

inline auto testMatrix2() {
  return Eigen::Matrix3d{
      {1.15,    0, 0},
      {   0, 0.75, 0},
      {   0,    0, 1}
  };
}

auto testMatPar() {
  double Emod = 1000;
  double nu   = 0.25;
  return YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
}

// Check if the implementation of using functions from Eigen library is correct with regard to explicit implementation
auto checkTwoPointMaterialTensor() {
  using StrainTags::deformationGradient;

  TestSuite t("Check Two-Point Material Tensor Implementation");
  spdlog::info("Testing " + t.name());

  auto matPar = testMatPar();
  auto mu     = convertLameConstants(matPar).toShearModulus();
  auto Lambda = convertLameConstants(matPar).toLamesFirstParameter();
  auto K      = convertLameConstants(matPar).toBulkModulus();

  std::array<double, 1> mu_og    = {mu};
  std::array<double, 1> alpha_og = {2.0};

  auto nh = Materials::NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

  constexpr double tol = 1e-15;

  auto F              = testMatrix1();
  auto pk2Stress      = nh.template stresses<deformationGradient, false>(F);
  auto materialTensor = nh.template tangentModuli<deformationGradient, false>(F);

  constexpr int dim                    = 3;
  const Eigen::Matrix<double, 3, 3> Id = Eigen::Matrix<double, 3, 3>::Identity();
  Eigen::TensorFixedSize<double, Eigen::Sizes<dim, dim, dim, dim>> A_exp;
  A_exp.setZero();
  for (const auto i : Dune::range(dim))
    for (const auto J : Dune::range(dim))
      for (const auto k : Dune::range(dim))
        for (const auto L : Dune::range(dim))
          for (const auto I : Dune::range(dim))
            for (const auto K : Dune::range(dim)) {
              A_exp(i, J, k, L) += materialTensor(I, J, K, L) * F(i, I) * F(k, K) + Id(i, k) * pk2Stress(J, L);
            }

  const auto A =
      transformTangentModuli<TangentModuliTags::Material, TangentModuliTags::TwoPoint>(materialTensor, pk2Stress, F);

  for (const auto i : Dune::range(dim))
    for (const auto J : Dune::range(dim))
      for (const auto k : Dune::range(dim))
        for (const auto L : Dune::range(dim)) {
          const std::string& indexName = "i = " + std::to_string(i) + ", J = " + std::to_string(J) +
                                         ", k = " + std::to_string(k) + ", L = " + std::to_string(L);
          checkScalars(t, A(i, J, k, L), A_exp(i, J, k, L), " Incorrect entry for two-point tensor at " + indexName,
                       tol);
        }

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  checkTwoPointMaterialTensor();

  return t.exit();
}
