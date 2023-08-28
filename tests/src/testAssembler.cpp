// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testHelpers.hh"

#include <vector>

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

using Dune::TestSuite;
#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/mechanics/nonLinearElastic.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

template <typename TestSuiteType, typename SparseType, typename DenseType, typename DOFSize>
void checkAssembledQuantities(TestSuiteType& t, SparseType& sType, DenseType& dType, DOFSize dofSize) {
  t.check(isApproxSame(sType, dType, 1e-15), "Dense==Sparse");
  t.check(sType.rows() == dofSize) << "DOFsCheck via rows: " << sType.rows() << "rows and " << dofSize << " DOFs";
  if (not(std::is_same_v<SparseType, Eigen::VectorXd>))
    t.check(sType.cols() == dofSize) << "DOFsCheck via columns: " << sType.cols() << "cols and " << dofSize << " DOFs";
}

auto SimpleAssemblersTest() {
  TestSuite t("SimpleAssemblersTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  for (int i = 0; i < 4; ++i) {
    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    auto basis        = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));
    auto matParameter = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});

    Ikarus::StVenantKirchhoff matSVK(matParameter);
    auto reducedMat = planeStress(matSVK, 1e-8);
    std::vector<Ikarus::NonLinearElastic<decltype(basis), decltype(reducedMat)>> fes;

    auto volumeLoad = []([[maybe_unused]] const auto& globalCoord, const auto& lamb) {
      Eigen::Vector2d fext;
      fext.setZero();
      fext[1] = 2 * lamb;
      fext[0] = lamb;
      return fext;
    };
    for (auto&& ge : elements(gridView))
      fes.emplace_back(basis, ge, reducedMat, volumeLoad);

    auto basisP = std::make_shared<const decltype(basis)>(basis);
    Ikarus::DirichletValues dirichletValues(basisP->flat());
    dirichletValues.fixDOFs([](auto& basis_, auto& dirichletFlags) {
      Dune::Functions::forEachBoundaryDOF(basis_, [&](auto&& indexGlobal) { dirichletFlags[indexGlobal] = true; });
    });

    Ikarus::SparseFlatAssembler sparseFlatAssembler(fes, dirichletValues);
    Ikarus::DenseFlatAssembler denseFlatAssembler(fes, dirichletValues);

    Eigen::VectorXd d(basis.flat().size());
    d.setRandom();
    double load = 0.0;

    Ikarus::FErequirements req = Ikarus::FErequirements()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, d)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, load)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness);

    auto& KRawDense = denseFlatAssembler.getRawMatrix(req);
    auto& KRaw      = sparseFlatAssembler.getRawMatrix(req);
    auto& RRawDense = denseFlatAssembler.getRawVector(req);
    auto& RRaw      = sparseFlatAssembler.getRawVector(req);
    checkAssembledQuantities(t, KRaw, KRawDense, 2 * gridView.size(2));
    checkAssembledQuantities(t, RRaw, RRawDense, 2 * gridView.size(2));

    auto& KDense = denseFlatAssembler.getMatrix(req);
    auto& K      = sparseFlatAssembler.getMatrix(req);
    auto& RDense = denseFlatAssembler.getVector(req);
    auto& R      = sparseFlatAssembler.getVector(req);
    checkAssembledQuantities(t, K, KDense, 2 * gridView.size(2));
    checkAssembledQuantities(t, R, RDense, 2 * gridView.size(2));

    const auto fixedDOFs    = dirichletValues.fixedDOFsize();
    const int boundaryNodes = (elementsPerDirection[0] * Dune::power(2, i) + 1) * 2
                              + (elementsPerDirection[1] * Dune::power(2, i) + 1) * 2 - 4;
    t.check(2 * boundaryNodes == fixedDOFs);

    auto& KRedDense = denseFlatAssembler.getReducedMatrix(req);
    auto& KRed      = sparseFlatAssembler.getReducedMatrix(req);
    auto& RRedDense = denseFlatAssembler.getReducedVector(req);
    auto& RRed      = sparseFlatAssembler.getReducedVector(req);
    checkAssembledQuantities(t, KRed, KRedDense, 2 * gridView.size(2) - fixedDOFs);
    checkAssembledQuantities(t, RRed, RRedDense, 2 * gridView.size(2) - fixedDOFs);

    grid->globalRefine(1);
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(SimpleAssemblersTest());
  return t.exit();
}
