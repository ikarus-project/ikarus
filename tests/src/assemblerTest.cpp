// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>

#include "testHelpers.hh"

#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
using Dune::TestSuite;
#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/mechanics/nonLinearElasticityFE.hh>

auto SimpleAssemblersTest() {
  TestSuite t("SimpleAssemblersTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  for (int i = 0; i < 4; ++i) {
    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

    std::vector<Ikarus::NonLinearElasticityFE<decltype(basis)>> fes;
    const double Emodul = 1000;
    auto volumeLoad
        = []([[maybe_unused]] const auto& globalCoord, const auto& lamb) {  // FIXME makeAnalytic globa function
            Eigen::Vector2d fext;
            fext.setZero();
            fext[1] = 2 * lamb;
            fext[0] = lamb;
            return fext;
          };
    for (auto&& ge : elements(gridView))
      fes.emplace_back(basis, ge, Emodul, 0.3, nullptr, nullptr, volumeLoad);

    std::vector<bool> dirichFlags(basis.size(), false);

    Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

    Ikarus::SparseFlatAssembler sparseFlatAssembler(basis, fes, dirichFlags);
    Ikarus::DenseFlatAssembler denseFlatAssembler(basis, fes, dirichFlags);

    Eigen::VectorXd d(basis.size());
    d.setRandom();
    Ikarus::FErequirements req = Ikarus::FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, d)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, 0)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                     .build();

    auto& Kdense = denseFlatAssembler.getMatrix(req);
    auto& K      = sparseFlatAssembler.getMatrix(req);

    const auto fixedDofs = std::ranges::count(dirichFlags, true);
    t.check(isApproxSame(K, Kdense, 1e-15), "Dense==Sparse");
    t.check(K.rows() == 2 * gridView.size(2), "DofsCheck");
    t.check(K.cols() == 2 * gridView.size(2), "DofsCheck");
    const int boundaryNodes = (elementsPerDirection[0] * Dune::power(2, i) + 1) * 2
                              + (elementsPerDirection[1] * Dune::power(2, i) + 1) * 2 - 4;
    t.check(2 * boundaryNodes == fixedDofs);

    auto& KdenseRed = denseFlatAssembler.getReducedMatrix(req);
    auto& KRed      = sparseFlatAssembler.getReducedMatrix(req);

    t.check(isApproxSame(KRed, KdenseRed, 1e-15), "DenseRed==SparseRed");
    t.check(KRed.rows() == 2 * gridView.size(2) - fixedDofs, "DofsCheckRed");
    t.check(KRed.cols() == 2 * gridView.size(2) - fixedDofs, "DofsCheckRed");

    grid->globalRefine(1);
  }
  return t;
}

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  t.subTest(SimpleAssemblersTest());
  return t.exit();
}
