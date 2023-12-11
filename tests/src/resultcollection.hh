// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "testcommon.hh"

#include <ikarus/finiteelements/ferequirements.hh>

template <typename FiniteElement>
auto getVertexPositions(FiniteElement& fe) {
  constexpr int dim            = FiniteElement::Traits::mydim;
  const auto& element          = fe.gridElement();
  const auto& referenceElement = Dune::referenceElement<double, dim>(element.type());
  const int numberOfVertices   = referenceElement.size(dim);

  std::vector<typename FiniteElement::GridElementEntityType::Geometry::LocalCoordinate> positions;
  for (auto i : std::views::iota(0, numberOfVertices))
    positions.push_back(referenceElement.position(i, dim));

  return positions;
}

inline auto linearStressResultsOfSquare = [](auto& nonLinearOperator, auto& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  if constexpr (requires { fe.setEASType(4); }) {
    fe.setEASType(4);
    expectedStress << 1214.28571429, 1214.28571429, 384.61538462, 1214.28571429, 214.28571429, 384.61538462,
        214.28571429, 1214.28571429, 384.61538462, 214.28571429, 214.28571429, 384.61538462;
  } else {
    expectedStress << 1428.57142857, 1428.57142857, 769.23076923, 1098.90109890, 329.67032967, 384.61538462,
        329.67032967, 1098.90109890, 384.61538462, 0, 0, 0;
  }

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 1, 1, 1, 1, 1, 1;

  using namespace Ikarus;
  auto resultRequirements = ResultRequirements()
                                .insertGlobalSolution(FESolutions::displacement, displacement)
                                .addResultRequest(ResultType::linearStress);

  return std::make_tuple(resultRequirements, expectedStress, getVertexPositions(fe));
};

inline auto linearStressResultsOfCube = [](auto& nonLinearOperator, auto& fe) {
  constexpr int vertices   = 8;
  constexpr int quantities = 6;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  expectedStress << 576.92307692, 1346.15384615, 576.92307692, 384.61538462, 0, 384.61538462, 0, 0, 0, 0, 0, 0,
      -1346.15384615, 192.30769231, -1346.15384615, 0, -769.23076923, 0, -1346.15384615, -576.92307692, -576.92307692,
      0, -384.61538462, -384.61538462, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -576.92307692, -576.92307692, -1346.15384615,
      -384.61538462, -384.61538462, 0, 0, 0, 0, 0, 0, 0;

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  using namespace Ikarus;
  auto resultRequirements = ResultRequirements()
                                .insertGlobalSolution(FESolutions::displacement, displacement)
                                .addResultRequest(ResultType::linearStress);

  return std::make_tuple(resultRequirements, expectedStress, getVertexPositions(fe));
};

inline auto linearStressResultsOfTriangle = [](auto& nonLinearOperator, auto& fe) {
  constexpr int vertices   = 3;
  constexpr int quantities = 3;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  expectedStress.rowwise() = Eigen::Matrix<double, 1, quantities>{2197.80219780, 659.34065934, 384.61538462};

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 2, 0, 1, 0;

  using namespace Ikarus;
  auto resultRequirements = ResultRequirements()
                                .insertGlobalSolution(FESolutions::displacement, displacement)
                                .addResultRequest(ResultType::linearStress);

  return std::make_tuple(resultRequirements, expectedStress, getVertexPositions(fe));
};

inline auto linearStressResultsOfTetrahedron = [](auto& nonLinearOperator, auto& fe) {
  constexpr int vertices   = 4;
  constexpr int quantities = 6;

  Eigen::Matrix<double, vertices, quantities> expectedStress;
  expectedStress.rowwise()
      = Eigen::Matrix<double, 1, quantities>{576.92307692, 1346.15384615, 576.92307692, 0, 384.61538462, 769.23076923};

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0;

  using namespace Ikarus;
  auto resultRequirements = ResultRequirements()
                                .insertGlobalSolution(FESolutions::displacement, displacement)
                                .addResultRequest(ResultType::linearStress);

  return std::make_tuple(resultRequirements, expectedStress, getVertexPositions(fe));
};