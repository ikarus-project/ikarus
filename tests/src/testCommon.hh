// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <vector>

#include <dune/alugrid/grid.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/nurbsgrid.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>

#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feBases/powerBasisFE.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/utils/duneUtilities.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/functionSanityChecks.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Grids {
  struct Yasp {};
  struct Alu {};
  struct Iga {};
}  // namespace Grids

template <typename GridType>
auto createGrid([[maybe_unused]] int elex = 10, [[maybe_unused]] int eley = 10) {
  //  //  /// ALUGrid Example
  if constexpr (std::is_same_v<GridType, Grids::Alu>) {
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
    auto grid  = Dune::GmshReader<Grid>::read("testFiles/unstructuredTrianglesfine.msh", false);
    grid->globalRefine(0);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::Yasp>) {
    using Grid     = Dune::YaspGrid<2>;
    const double L = 1;
    const double h = 1;

    Dune::FieldVector<double, 2> bbox       = {L, h};
    std::array<int, 2> elementsPerDirection = {elex, eley};
    auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::Iga>) {
    constexpr auto dimworld        = 2;
    const std::array<int, 2> order = {2, 2};

    const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

    using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

    const std::vector<std::vector<ControlPoint>> controlPoints
        = {{{.p = {0, 0}, .w = 5}, {.p = {0.5, 0}, .w = 1}, {.p = {1, 0}, .w = 1}},
           {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 10}, {.p = {1, 0.5}, .w = 1}},
           {{.p = {0, 1}, .w = 1}, {.p = {0.5, 1}, .w = 1}, {.p = {1, 1}, .w = 1}}};

    std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

    auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
    using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

    Dune::IGA::NURBSPatchData<2, dimworld> patchData;
    patchData.knotSpans     = knotSpans;
    patchData.degree        = order;
    patchData.controlPoints = controlNet;
    auto grid               = std::make_shared<Grid>(patchData);
    grid->globalRefine(1);
    return grid;
  }
}

template <int size>
struct CornerFactory {
  static void construct(std::vector<Dune::FieldVector<double, size>>& values, const int corners = 10) {
    values.resize(corners);
    std::generate(values.begin(), values.end(),
                  []() { return Ikarus::createRandomVector<Dune::FieldVector<double, size>>(); });
  }
};

enum class CornerDistortionFlag { unDistorted, fixedDistorted, randomlyDistorted };

// Corner factory for element with codim==0, e.g. no surfaces in 3D
template <int gridDim>
struct ValidCornerFactory {
  static void construct(std::vector<Dune::FieldVector<double, gridDim>>& values, const Dune::GeometryType& type,
                        const CornerDistortionFlag& distortionFlag) {
    const auto& refElement = Dune::ReferenceElements<double, gridDim>::general(type);

    const auto numberOfVertices = refElement.size(gridDim);

    values.resize(numberOfVertices);
    for (int i = 0; i < numberOfVertices; ++i)
      values[i] = refElement.position(i, gridDim);

    /// Perturbation of the corner values
    if (distortionFlag == CornerDistortionFlag::unDistorted)
      return;
    else if (distortionFlag == CornerDistortionFlag::fixedDistorted) {
      if (not((gridDim == 2) and (numberOfVertices == 4)))
        DUNE_THROW(Dune::NotImplemented, "Fixed distortion is only implemented for a 2D 4-node element (Q1)");
      std::vector<Dune::FieldVector<double, gridDim>> randomnessOnNodes;
      randomnessOnNodes.push_back({-0.2, -0.05});
      randomnessOnNodes.push_back({-0.15, 0.05});
      randomnessOnNodes.push_back({0.15, 0.15});
      randomnessOnNodes.push_back({-0.05, -0.1});
      for (long i = 0; i < 4; ++i)
        values[i] += randomnessOnNodes[i];
    } else if (distortionFlag == CornerDistortionFlag::randomlyDistorted) {
      std::transform(values.begin(), values.end(), values.begin(), [](const auto& vec) {
        return vec + Ikarus::createRandomVector<Dune::FieldVector<double, gridDim>>(-0.2, 0.2);
      });
    } else
      DUNE_THROW(Dune::IOError, "Incorrect CornerDistortionFlag entered");
  }
};

template <int gridDim>
auto createUGGridFromCorners(const CornerDistortionFlag& distortionFlag) {
  using Grid = Dune::UGGrid<gridDim>;

  std::vector<Dune::FieldVector<double, gridDim>> corners;

  const int numberOfVertices = Dune::power(2, gridDim);
  ValidCornerFactory<gridDim>::construct(corners, Dune::GeometryTypes::cube(gridDim), distortionFlag);

  std::vector<unsigned int> vertexArrangment;
  vertexArrangment.resize(numberOfVertices);
  std::iota(vertexArrangment.begin(), vertexArrangment.end(), 0);

  Dune::GridFactory<Grid> gridFactory;
  for (auto& corner : corners) {
    gridFactory.insertVertex(corner);
  }
  gridFactory.insertElement(Dune::GeometryTypes::cube(gridDim), vertexArrangment);

  std::unique_ptr<Grid> grid = gridFactory.createGrid();
  return grid;
}

template <typename Ele>
struct ElementTest {};

template <typename NonLinearOperator>
[[nodiscard]] auto checkGradientOfElement(NonLinearOperator& nonLinearOperator,
                                          const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check gradient");
  t.check(checkGradient(nonLinearOperator, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "The gradient of calculateVector is not the gradient of calculateScalar." << messageIfFailed;
  return t;
}

template <typename NonLinearOperator>
[[nodiscard]] auto checkHessianOfElement(NonLinearOperator& nonLinearOperator,
                                         const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check Hessian");
  t.check(checkHessian(nonLinearOperator, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "The Hessian of calculateMatrix is not the Hessian of calculateScalar. " << messageIfFailed;
  return t;
}

template <typename NonLinearOperator>
[[nodiscard]] auto checkJacobianOfElement(NonLinearOperator& nonLinearOperator,
                                          const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check Jacobian");
  t.check(checkJacobian(nonLinearOperator, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "The Jacobian of calculateMatrix is not the Jacobian of calculateVector." << messageIfFailed;
  return t;
}

template <typename NonLinearOperator, typename FiniteElement>
[[nodiscard]] auto checkCauchyStressOf2DElement(NonLinearOperator& nonLinearOperator, FiniteElement& fe,
                                                const std::string& messageIfFailed = "") {
  static_assert(FiniteElement::Traits::mydim == 2,
                "The test to check Cauchy stress is only supported for the two-dimensional case");
  assert(fe.localView().element().type() == Dune::GeometryTypes::quadrilateral
         && "The test to check Cauchy stress is only supported for quadrilaterals");
  static_assert(std::remove_cvref_t<decltype(fe.localView().tree())>::degree() == 2,
                "The test to check Cauchy stress is only supported with the powerBasis being two-dimensional");
  static_assert(std::is_same_v<typename std::remove_cvref_t<decltype(fe.localView().tree().child(0))>,
                               Dune::Functions::LagrangeNode<
                                   std::remove_cvref_t<decltype(fe.localView().globalBasis().gridView())>, 1, double>>,
                "The test to check Cauchy stress is only supported for a linear Lagrange basis");

  Dune::TestSuite t("Cauchy Stress check of 2D solid element (4-node quadrilateral)");
  using namespace Ikarus;
  using namespace Dune::Indices;
  using namespace Dune::Functions::BasisFactory;
  const int gridDim = 2;

  auto& displacement = nonLinearOperator.firstParameter();
  displacement << 0, 0, 1, 1, 1, 1, 1, 1;
  std::array<Dune::FieldVector<double, 3>, 4> expectedStress = {{{1428.5714287000, 1428.5714287000, 769.2307693000},
                                                                 {1098.9010990000, 329.6703297000, 384.6153846500},
                                                                 {329.6703297000, 1098.9010990000, 384.6153846500},
                                                                 {0, 0, 0}}};

  /// Modifications for the EAS method
  /// The test is only to check the computation of Cauchy stress in the linear elastic case. This includes only the
  /// addition of C*M*alpha to the stress from pure displacement formulation. Hence it is sufficient to test only for
  /// one of the EAS variant as the other M matrices are already tested.
  if constexpr (requires { fe.setEASType(4); }) {
    fe.setEASType(4);
    expectedStress = {{{1214.2857143755, 1214.2857143755, 384.6153846800},
                       {1214.2857144055, 214.2857142945, 384.6153847400},
                       {214.2857142945, 1214.2857144055, 384.6153845600},
                       {214.2857143245, 214.2857143245, 384.6153846200}}};
  }

  auto resultRequirements = Ikarus::ResultRequirements<>()
                                .insertGlobalSolution(Ikarus::FESolutions::displacement, displacement)
                                .addResultRequest(ResultType::linearStress);

  ResultTypeMap<double> result;
  auto gridView        = fe.localView().globalBasis().gridView();
  auto scalarBasis     = makeConstSharedBasis(gridView, lagrangeDG<1>());
  auto localScalarView = scalarBasis->localView();
  std::vector<Dune::FieldVector<double, 3>> stressVector(scalarBasis->size());
  auto ele = elements(gridView).begin();

  localScalarView.bind(*ele);
  const auto& fe2              = localScalarView.tree().finiteElement();
  const auto& referenceElement = Dune::ReferenceElements<double, gridDim>::general(ele->type());
  for (auto c = 0UL; c < fe2.size(); ++c) {
    const auto fineKey                        = fe2.localCoefficients().localKey(c);
    const auto nodalPositionInChildCoordinate = referenceElement.position(fineKey.subEntity(), fineKey.codim());
    fe.calculateAt(resultRequirements, nodalPositionInChildCoordinate, result);
    Eigen::Vector3d computedResult = result.getResult(ResultType::linearStress);
    const auto nodeIndex           = localScalarView.index(localScalarView.tree().localIndex(c))[0];
    stressVector[nodeIndex]        = toDune(computedResult);
    for (auto voigtIndex = 0UL; voigtIndex < 3; ++voigtIndex) {
      const auto FEStressComponent       = stressVector[nodeIndex][voigtIndex];
      const auto expectedStressComponent = expectedStress[c][voigtIndex];
      const bool isStressCorrect
          = Dune::FloatCmp::eq<double, Dune::FloatCmp::CmpStyle::absolute>(FEStressComponent, expectedStressComponent);
      t.check(isStressCorrect) << "Stress component " << voigtIndex << " at node " << c << " is not correctly computed"
                               << messageIfFailed;
    }
  }
  return t;
}

template <typename NonLinearOperator, typename FiniteElement,
          typename FERequirementType = FiniteElement::FERequirementType>
[[nodiscard]] auto checkFEByAutoDiff(NonLinearOperator&, FiniteElement& fe, FERequirementType req,
                                     const std::string& messageIfFailed = "") {
  Dune::TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation");
  auto& basis           = fe.localView().globalBasis();
  auto nDOF             = basis.size();
  using AutoDiffBasedFE = Ikarus::AutoDiffFE<FiniteElement, FERequirementType, false, true>;
  AutoDiffBasedFE feAutoDiff(fe);

  const double tol = 1e-10;

  Eigen::MatrixXd K, KAutoDiff;
  K.setZero(nDOF, nDOF);
  KAutoDiff.setZero(nDOF, nDOF);

  fe.calculateMatrix(req, K);
  feAutoDiff.calculateMatrix(req, KAutoDiff);

  if constexpr (requires { feAutoDiff.getFE().getNumberOfEASParameters(); }) {
    t.check(fe.getNumberOfEASParameters() == feAutoDiff.getFE().getNumberOfEASParameters())
        << "Number of EAS parameters for FE(" << fe.getNumberOfEASParameters()
        << ") and number of EAS parameters for AutodiffFE(" << feAutoDiff.getFE().getNumberOfEASParameters()
        << ") are not equal";
  }

  t.check(K.isApprox(KAutoDiff, tol),
          "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
          "automatic differentiation")
      << messageIfFailed << "\nKAutoDiff:\n"
      << KAutoDiff << "\nK:\n"
      << K;

  try {
    Eigen::VectorXd R, RAutoDiff;
    R.setZero(nDOF);
    RAutoDiff.setZero(nDOF);

    fe.calculateVector(req, R);
    feAutoDiff.calculateVector(req, RAutoDiff);
    t.check(R.isApprox(RAutoDiff, tol),
            "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
            "automatic differentiation")
        << messageIfFailed << "\nRAutoDiff:\n"
        << RAutoDiff << "\nR:\n"
        << R;
  } catch (const Dune::NotImplemented&) {
    /// calculateVector not checked since the element throws a Dune::NotImplemented
  }

  try {
    auto energy         = fe.calculateScalar(req);
    auto energyAutoDiff = feAutoDiff.calculateScalar(req);
    t.check(Dune::FloatCmp::eq(energy, energyAutoDiff, tol),
            "Mismatch between the energies obtained from explicit implementation and the one based on "
            "automatic differentiation")
        << messageIfFailed;
  } catch (const Dune::NotImplemented&) {
  }

  return t;
}
