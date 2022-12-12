// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "common.hh"

#include <dune/common/bitsetvector.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/grid/uggrid.hh>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>

/** These tests test your element on some gridElement with some basis
 *
 * @tparam FEElementTemplate The element as template template parameter. The template needs to be the globalBasis
 * @tparam gridDim The dimension of the grid element the finite element should be tested
 * @tparam PreBasis The preBasis you want to test the element with
 * @tparam F A variadic number of the test functor you want to be checked, they need to accept a non-linear operator and
 * the finite element
 */
template <template <typename> typename FEElementTemplate, int gridDim, typename PreBasis, typename... F>
auto testFEElement(const PreBasis& preBasis, const std::string& elementName, F&&... f) {
  TestSuite t(std::string("testFEElement ") + elementName + " on grid element with dimension" + std::to_string(gridDim)
              + ".");

  auto fTuple = std::forward_as_tuple(f...);

  using Grid = Dune::UGGrid<gridDim>;

  std::vector<Dune::FieldVector<double, gridDim>> corners;

  const int numberOfVertices = Dune::power(2, gridDim);
  ValidCornerFactory<gridDim>::construct(corners, Dune::GeometryTypes::cube(gridDim));

  std::vector<unsigned int> vertexArrangment;
  vertexArrangment.resize(numberOfVertices);
  std::iota(vertexArrangment.begin(), vertexArrangment.end(), 0);

  Dune::GridFactory<Grid> gridFactory;
  for (auto& corner : corners) {
    gridFactory.insertVertex(corner);
  }
  gridFactory.insertElement(Dune::GeometryTypes::cube(gridDim), vertexArrangment);

  std::unique_ptr<Grid> grid = gridFactory.createGrid();

  auto gridView = grid->leafGridView();
  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, preBasis);

  auto localView = basis.localView();

  auto volumeLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    VectorType fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };

  auto neumannBoundaryLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    VectorType fext;
    fext.setZero();
    fext[1] = lamb / 40;
    fext[0] = 0;
    return fext;
  };

  // We artificially apply a Neumann load on the complete boundary
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  const double youngsModulus = 1000;
  const double poissonsRatio = 0.3;
  auto element               = gridView.template begin<0>();

  using FEElementType = FEElementTemplate<decltype(basis)>;
  std::vector<FEElementType> fes;
  fes.emplace_back(basis, *element, youngsModulus, poissonsRatio, &volumeLoad, &neumannBoundary, &neumannBoundaryLoad);
  auto basisP = std::make_shared<const decltype(basis)>(basis);

  Ikarus::DirichletValues dirichletValues(basisP);
  auto& fe             = fes[0];
  auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

  typename FEElementType::FERequirementType::SolutionVectorType d;
  d.setRandom(basis.size());

  double lambda = 7.3;

  auto requirements = FErequirements()
                          .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                          .addAffordance(Ikarus::AffordanceCollections::elastoStatics);
  Eigen::VectorXd forces;
  Eigen::MatrixXd stiffnessmatrix;

  auto fvLambda = [&](auto&& d_) -> auto {
    forces.setZero(basis.localView().maxSize());
    requirements.insertGlobalSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getScalar(requirements);
  };

  auto dfvLambda = [&](auto&& d_) -> auto& {
    forces.setZero(basis.localView().maxSize());
    requirements.insertGlobalSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getVector(requirements);
  };
  auto ddfvLambda = [&](auto&& d_) -> auto& {
    requirements.insertGlobalSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getMatrix(requirements);
  };
  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(fvLambda, dfvLambda, ddfvLambda), parameter(d));

  // execute all passed functions
  nonLinOp.updateAll();
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<sizeof...(F)>()),
                        [&](auto i) { t.subTest(std::get<i.value>(fTuple)(nonLinOp, fe)); });

  // check if element has a test functor, if yes we execute it
  if constexpr (requires { ElementTest<FEElementType>::test(); }) {
    auto testFunctor = ElementTest<FEElementType>::test();
    t.subTest(testFunctor(nonLinOp, fe));
  }

  return t;
}

auto checkGradientFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe) { return checkGradientOfElement(nonLinOp); };
auto checkHessianFunctor  = [](auto& nonLinOp, [[maybe_unused]] auto& fe) { return checkHessianOfElement(nonLinOp); };
auto checkJacobianFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe) {
  auto subOperator = nonLinOp.template subOperator<1, 2>();
  return checkJacobianOfElement(subOperator);
};
