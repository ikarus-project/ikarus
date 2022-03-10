//
// Created by Alex on 21.04.2021.
//
#include <gmock/gmock.h>

#include "../../config.h"
#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>

#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/NonLinearElasticityFEwithBasis.h>

TEST(Assembler, SimpleAssemblersTest) {
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox = {4, 2};
  std::array<int, 2> eles           = {2, 1};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  for (int i = 0; i < 4; ++i) {
    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

    const auto& indexSet = gridView.indexSet();

    std::vector<Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(basis)>> fes;
    const double Emodul = 1000;
    auto volumeLoad = [](const auto& globalCoord, const auto& lamb)
    {Eigen::Vector2d fext;
      fext.setZero();
      fext[1] = 2 * lamb;
      fext[0] = lamb;
      return fext;
    };
    for (auto&& ge : elements(gridView))
      fes.emplace_back(basis, ge, Emodul, 0.3,volumeLoad);

    std::vector<bool> dirichFlags(basis.size(), false);

    Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

    Ikarus::SparseFlatAssembler sparseFlatAssembler(basis, fes, dirichFlags);
    Ikarus::DenseFlatAssembler denseFlatAssembler(basis, fes, dirichFlags);

    Ikarus::FErequirements<Eigen::VectorXd> req;
    Eigen::VectorXd d(basis.size());
    d.setRandom();
    req.sols.emplace_back(d);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, 0});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    auto& Kdense          = denseFlatAssembler.getMatrix(req);
    auto& K               = sparseFlatAssembler.getMatrix(req);

    const auto fixedDofs = std::ranges::count(dirichFlags, true);
    EXPECT_THAT(K, EigenApproxEqual(Kdense, 1e-15));
    EXPECT_THAT(K.rows(), 2 * gridView.size(2));
    EXPECT_THAT(K.cols(), 2 * gridView.size(2));
    const int boundaryNodes = (eles[0] * Dune::power(2, i) + 1) * 2 + (eles[1] * Dune::power(2, i) + 1) * 2 - 4;
    EXPECT_EQ(fixedDofs, 2 * boundaryNodes);

    auto& KdenseRed = denseFlatAssembler.getReducedMatrix(req);
    auto& KRed      = sparseFlatAssembler.getReducedMatrix(req);

    EXPECT_THAT(KRed, EigenApproxEqual(KdenseRed, 1e-15));
    EXPECT_THAT(KRed.rows(), 2 * gridView.size(2) - fixedDofs);
    EXPECT_THAT(KRed.cols(), 2 * gridView.size(2) - fixedDofs);

    grid->globalRefine(1);
  }
}
