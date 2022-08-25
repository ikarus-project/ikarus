

#include <config.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "testHelpers.hh"

#include <vector>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/mechanics/nonLinearElasticityFE.hh>

TEST_CASE("Assembler: SimpleAssemblersTest", "[assemblerTest.cpp]") {
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox = {4, 2};
  std::array<int, 2> eles           = {2, 1};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  for (int i = 0; i < 4; ++i) {
    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

    const auto& indexSet = gridView.indexSet();

    std::vector<Ikarus::NonLinearElasticityFE<decltype(basis)>> fes;
    const double Emodul = 1000;
    auto volumeLoad     = [](const auto& globalCoord, const auto& lamb) {  // FIXME makeAnalytic globa function
      Eigen::Vector2d fext;
      fext.setZero();
      fext[1] = 2 * lamb;
      fext[0] = lamb;
      return fext;
    };
    for (auto&& ge : elements(gridView))
      fes.emplace_back(basis, ge, Emodul, 0.3, volumeLoad);

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
    CHECK_THAT(K, EigenApproxEqual(Kdense, 1e-15));
    CHECK(K.rows() == 2 * gridView.size(2));
    CHECK(K.cols() == 2 * gridView.size(2));
    const int boundaryNodes = (eles[0] * Dune::power(2, i) + 1) * 2 + (eles[1] * Dune::power(2, i) + 1) * 2 - 4;
    CHECK(2 * boundaryNodes == fixedDofs);

    auto& KdenseRed = denseFlatAssembler.getReducedMatrix(req);
    auto& KRed      = sparseFlatAssembler.getReducedMatrix(req);

    CHECK_THAT(KRed, EigenApproxEqual(KdenseRed, 1e-15));
    CHECK(KRed.rows() == 2 * gridView.size(2) - fixedDofs);
    CHECK(KRed.cols() == 2 * gridView.size(2) - fixedDofs);

    grid->globalRefine(1);
  }
}
