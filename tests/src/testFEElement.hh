// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testCommon.hh"

#include <dune/common/bitsetvector.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/grid/uggrid.hh>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/io/resultFunction.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/utils/basis.hh>

/** These tests test your element on some gridElement with some basis
 *
 * @tparam FEElementTemplate The element as template template parameter. The template needs to be the globalBasis
 * @tparam gridDim The dimension of the grid element the finite element should be tested
 * @tparam PreBasis The preBasis you want to test the element with
 * @tparam F A variadic number of the test functor you want to be checked, they need to accept a non-linear operator and
 * the finite element
 */
template <template <typename> typename FEElementTemplate, int gridDim, typename PreBasis, typename... F>
auto testFEElement(const PreBasis& preBasis, const std::string& elementName, const CornerDistortionFlag& distortionFlag,
                   F&&... f) {
  Dune::TestSuite t(std::string("testFEElement ") + elementName + " on grid element with dimension "
                    + std::to_string(gridDim));

  auto fTuple = std::forward_as_tuple(f...);

  auto grid = createUGGridFromCorners<gridDim>(distortionFlag);

  auto gridView = grid->leafGridView();
  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis     = Ikarus::makeBasis(gridView, preBasis);
  auto flatBasis = basis.flat();

  auto localView = flatBasis.localView();

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
  fes.emplace_back(basis, *element, youngsModulus, poissonsRatio, volumeLoad, &neumannBoundary, neumannBoundaryLoad);
  auto basisP = std::make_shared<const decltype(basis)>(basis);

  Ikarus::DirichletValues dirichletValues(basisP->flat());
  auto& fe             = fes[0];
  auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

  static_assert(std::is_reference_v<typename decltype(sparseAssembler)::FEContainerType>);

  typename FEElementType::FERequirementType::SolutionVectorTypeRaw d;
  d.setRandom(basis.flat().size());

  double lambda = 7.3;

  auto requirements = FErequirements()
                          .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                          .addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto fvLambda = [&](auto&& d_) -> auto {
    requirements.insertGlobalSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getScalar(requirements);
  };

  auto dfvLambda = [&](auto&& d_) -> auto& {
    requirements.insertGlobalSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getVector(requirements);
  };
  auto ddfvLambda = [&](auto&& d_) -> auto& {
    requirements.insertGlobalSolution(Ikarus::FESolutions::displacement, d_);
    return sparseAssembler.getMatrix(requirements);
  };
  auto nonLinOp = Ikarus::NonLinearOperator(functions(fvLambda, dfvLambda, ddfvLambda), parameter(d));

  // execute all passed functions
  nonLinOp.updateAll();
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<sizeof...(F)>()),
                        [&](auto i) { t.subTest(std::get<i.value>(fTuple)(nonLinOp, fe, requirements)); });

  // check if element has a test functor, if yes we execute it
  if constexpr (requires { ElementTest<FEElementType>::test(); }) {
    auto testFunctor = ElementTest<FEElementType>::test();
    t.subTest(testFunctor(nonLinOp, fe, requirements));
  } else
    spdlog::info("No element test functor found for {}", Dune::className<FEElementType>());

  // Trying to instantiate the Result Evaluator
  if constexpr (gridDim == 2) {
    auto resReq         = Ikarus::ResultRequirements(requirements).addResultRequest(ResultType::PK2Stress);
    auto resultFunction = std::make_shared<ResultFunction<FEElementType>>(&fes, resReq);
  }

  return t;
}

auto checkGradientFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe, [[maybe_unused]] auto& req) {
  return checkGradientOfElement(nonLinOp);
};
auto checkHessianFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe, [[maybe_unused]] auto& req) {
  return checkHessianOfElement(nonLinOp);
};
auto checkJacobianFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe, [[maybe_unused]] auto& req) {
  auto subOperator = nonLinOp.template subOperator<1, 2>();
  return checkJacobianOfElement(subOperator);
};
auto checkCauchyStressFunctor
    = [](auto& nonLinOp, auto& fe, [[maybe_unused]] auto& req) { return checkCauchyStressOf2DElement(nonLinOp, fe); };
auto checkFEByAutoDiffFunctor
    = [](auto& nonLinOp, auto& fe, auto& req) { return checkFEByAutoDiff(nonLinOp, fe, req); };
