// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <ikarus/finiteelements/mechanics/materials/strainconversions.hh>
#include <ikarus/finiteelements/mechanics/materials/stressconversions.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;
using namespace Ikarus;

namespace Testing {
constexpr auto testDeformationGradient() {
  return Eigen::Matrix3d{
      {1.0, 0.0, 0.5},
      {0.0, 0.2, 0.5},
      {0.2, 0.5, 1.0}
  };
}

constexpr auto testStress() {
  return Eigen::Matrix3d{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };
}

auto testCauchyGreen() {
  auto F = Testing::testDeformationGradient();
  return (F.transpose() * F).eval();
}
} // namespace Testing

auto testStressConversions() {
  TestSuite t("Test Stress Conversions");
  const double epsilon = 1e-14;

  auto S = Testing::testStress();
  auto F = Testing::testDeformationGradient();

  auto stressTagRange = Dune::Hybrid::integralRange(std::integral_constant<int, 2>(), std::integral_constant<int, 5>());
  Dune::Hybrid::forEach(stressTagRange, [&](auto i) {
    constexpr StressTags tag1 = StressTags(int(i));
    Dune::Hybrid::forEach(stressTagRange, [&](auto j) {
      constexpr StressTags tag2         = StressTags(int(j));
      Eigen::Matrix3d convertedStress   = transformStress<tag1, tag2>(S, F);
      Eigen::Matrix3d reconvertedStress = transformStress<tag2, tag1>(convertedStress, F);
      checkApproxMatrices(
          t, reconvertedStress, S,
          "Re-converted stress not equal to initial stress (" + toString(tag1) + " & " + toString(tag2) + ")", epsilon);

      Dune::Hybrid::forEach(stressTagRange, [&](auto k) {
        constexpr StressTags tag3               = StressTags(int(k));
        Eigen::Matrix3d convertedAgainStress    = transformStress<tag2, tag3>(convertedStress, F);
        Eigen::Matrix3d directlyConvertedStress = transformStress<tag1, tag3>(S, F);

        checkApproxMatrices(t, convertedAgainStress, directlyConvertedStress,
                            "Converted stress with two steps does not equal direct conversion (" + toString(tag1) +
                                " & " + toString(tag2) + " & " + toString(tag3) + ")",
                            epsilon);
      });
    });
  });

  return t;
}

auto testStrainConversions() {
  TestSuite t("Test Strain Conversions");
  const double epsilon = 1e-14;

  Eigen::Matrix3d C = Testing::testCauchyGreen();

  auto strainRange = Dune::Hybrid::integralRange(std::integral_constant<int, 3>(), std::integral_constant<int, 5>());
  Dune::Hybrid::forEach(strainRange, [&](auto i) {
    constexpr StrainTags tag1 = StrainTags(int(i));
    Dune::Hybrid::forEach(strainRange, [&](auto j) {
      constexpr StrainTags tag2         = StrainTags(int(j));
      Eigen::Matrix3d convertedStrain   = transformStrain<tag1, tag2>(C);
      Eigen::Matrix3d reconvertedStrain = transformStrain<tag2, tag1>(convertedStrain);
      checkApproxMatrices(
          t, reconvertedStrain, C,
          "Re-converted strain not equal to initial strain (" + toString(tag1) + " & " + toString(tag2) + ")", epsilon);

      Dune::Hybrid::forEach(strainRange, [&](auto k) {
        constexpr StrainTags tag3               = StrainTags(int(k));
        Eigen::Matrix3d convertedAgainStrain    = transformStrain<tag2, tag3>(convertedStrain);
        Eigen::Matrix3d directlyConvertedStrain = transformStrain<tag1, tag3>(C);

        checkApproxMatrices(t, convertedAgainStrain, directlyConvertedStrain,
                            "Converted strain with two steps does not equal direct conversion (" + toString(tag1) +
                                " & " + toString(tag2) + " & " + toString(tag3) + ")",
                            epsilon);
      });
    });
  });

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testStressConversions());
  t.subTest(testStrainConversions());

  return t.exit();
}
