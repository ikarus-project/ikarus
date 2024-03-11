// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testcommon.hh"

#include <dune/common/bitsetvector.hh>
#include <dune/fufem/boundarypatch.hh>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/mechanics/loads.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/nonlinearoperator.hh>

/** These tests test your element on some gridElement with some basis
 *
 * \tparam PREFunc Type of the functor returning a pre finite element.
 * \tparam ReferenceElement The reference element, the grid has to be constructed from
 * \tparam PreBasis The preBasis you want to test the element with
 * \tparam F A variadic number of the test functor you want to be checked, they need to accept a non-linear operator and
 * the finite element
 */
template <typename PREFunc, typename PreBasis, typename ReferenceElement, typename Skills, typename... F>
auto testFEElement(const PreBasis& preBasis, const std::string& elementName, const CornerDistortionFlag& distortionFlag,
                   const ReferenceElement& refElement, PREFunc&& preFunc, Skills&& additionalSkills, F&&... f) {
  constexpr int gridDim = ReferenceElement::dimension;

  Dune::TestSuite t(std::string("testFEElement ") + elementName + " on grid element with dimension " +
                    std::to_string(gridDim));

  auto fTuple = std::forward_as_tuple(f...);

  auto grid     = createUGGridFromCorners<gridDim>(distortionFlag, refElement.type());
  auto gridView = grid->leafGridView();
  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis     = Ikarus::makeBasis(gridView, preBasis);
  auto flatBasis = basis.flat();

  auto localView = flatBasis.localView();

  auto vL = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    Eigen::Vector<typename VectorType::field_type, VectorType::dimension> fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };

  auto neumannBl = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    Eigen::Vector<typename VectorType::field_type, VectorType::dimension> fext;
    fext.setZero();
    fext[1] = lamb / 40;
    fext[0] = 0;
    return fext;
  };

  // We artificially apply a Neumann load on the complete boundary
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  auto element = gridView.template begin<0>();

  static constexpr int worldDim = std::remove_reference_t<decltype(*grid)>::dimensionworld;
  auto skillsMerged             = merge(skills(preFunc({.emodul = 1000, .nu = 0.3}), volumeLoad<worldDim>(vL),
                                               neumannBoundaryLoad(&neumannBoundary, neumannBl)),
                                        std::move(additionalSkills));
  using FEType                  = decltype(Ikarus::makeFE(basis, std::move(skillsMerged)));
  std::vector<FEType> fes;
  fes.emplace_back(Ikarus::makeFE(basis, std::move(skillsMerged)));
  fes[0].bind(*element);

  Ikarus::DirichletValues dirichletValues(basis.flat());
  auto& fe             = fes[0];
  auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

  static_assert(std::is_reference_v<typename decltype(sparseAssembler)::FEContainerType>);

  typename FEType::FERequirementType::SolutionVectorTypeRaw d;
  d.setRandom(basis.flat().size());

  double lambda = 7.3;

  auto requirements = FERequirements()
                          .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                          .addAffordance(Ikarus::AffordanceCollections::elastoStatics);
  bool successlinearStressCalculateAt = false;
  Dune::Hybrid::forEach(typename FEType::SupportedResultTypes(), [&, gridDim]<typename RT>(RT i) {
    static_assert(requires {
      fes[0].template calculateAt<RT::template Rebind>(requirements, Dune::FieldVector<double, gridDim>());
    });

    if ("linearStress" == toString(i))
      successlinearStressCalculateAt = true;
  });
  t.check(successlinearStressCalculateAt)
      << "linearStress call for calculateAt didn't work." << "\nThe supported types are "
      << Dune::className<typename FEType::SupportedResultTypes>() << "\n"
      << "The element is \n"
      << Dune::className<FEType>();

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
  if constexpr (requires { ElementTest<FEType>::test(); }) {
    auto testFunctor = ElementTest<FEType>::test();
    t.subTest(testFunctor(nonLinOp, fe, requirements));
  } else
    spdlog::info("No element test functor found for {}", Dune::className<FEType>());

  return t;
}

inline auto checkGradientFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe, [[maybe_unused]] auto& req) {
  return checkGradientOfElement(nonLinOp);
};
inline auto checkHessianFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe, [[maybe_unused]] auto& req) {
  return checkHessianOfElement(nonLinOp);
};
inline auto checkJacobianFunctor = [](auto& nonLinOp, [[maybe_unused]] auto& fe, [[maybe_unused]] auto& req) {
  auto subOperator = nonLinOp.template subOperator<1, 2>();
  return checkJacobianOfElement(subOperator);
};

template <template <typename, int, int> class RT, typename ResultEvaluator = Ikarus::Impl::DefaultUserFunction>
auto checkResultFunctionFunctorFactory(const auto& resultCollectionFunction) {
  return [&](auto& nonLinOp, auto& fe, [[maybe_unused]] auto& req) {
    auto [feRequirements, expectedStress, positions] = resultCollectionFunction(nonLinOp, fe);
    return checkResultFunction<RT, ResultEvaluator>(nonLinOp, fe, feRequirements, expectedStress, positions);
  };
}

inline auto checkFEByAutoDiffFunctor = [](auto& nonLinOp, auto& fe, auto& req) {
  return checkFEByAutoDiff(nonLinOp, fe, req);
};

template <template <typename, int, int> class RT, bool voigt = true>
auto checkCalculateAtFunctorFactory(const auto& resultCollectionFunction) {
  return [&](auto& nonLinOp, auto& fe, [[maybe_unused]] auto& req) {
    auto [feRequirements, expectedStress, positions] = resultCollectionFunction(nonLinOp, fe);
    if constexpr (voigt)
      return checkCalculateAt<RT>(nonLinOp, fe, feRequirements, expectedStress, positions);
    else
      return checkCalculateAt<RT, voigt>(nonLinOp, fe, feRequirements, stressResultsToMatrix(expectedStress),
                                         positions);
  };
}
