// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

using Dune::TestSuite;

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/svk.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

template <typename TestSuiteType, typename SparseType, typename DenseType>
void checkAssembledQuantities(TestSuiteType& t, SparseType& sType, DenseType& dType, std::size_t dofSize) {
  t.check(isApproxSame(sType, dType, 1e-15), "Dense==Sparse");
  t.check(sType.rows() == static_cast<Eigen::Index>(dofSize))
      << "DOFsCheck via rows: " << sType.rows() << "rows and " << dofSize << " DOFs";
  if (not(std::is_same_v<std::remove_cvref_t<SparseType>, Eigen::VectorXd>))
    t.check(sType.cols() == static_cast<Eigen::Index>(dofSize))
        << "DOFsCheck via columns: " << sType.cols() << " cols and " << dofSize << " DOFs";
}

template <typename PreBasis>
auto SimpleAssemblersTest(const PreBasis& preBasis) {
  TestSuite t("SimpleAssemblersTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  for (int ref = 0; ref < 4; ++ref) {
    auto gridView = grid->leafGridView();

    auto basis        = Ikarus::makeBasis(gridView, preBasis);
    auto matParameter = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});
    auto totalDOFs    = basis.flat().size();

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

    Ikarus::FERequirements req = Ikarus::FERequirements()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, d)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, load)
                                     .addAffordance(Ikarus::AffordanceCollections::elastoStatics);

    auto& KRawDense = denseFlatAssembler.getRawMatrix(req);
    auto& KRaw      = sparseFlatAssembler.getRawMatrix(req);
    auto& RRawDense = denseFlatAssembler.getRawVector(req);
    auto& RRaw      = sparseFlatAssembler.getRawVector(req);
    checkAssembledQuantities(t, KRaw, KRawDense, totalDOFs);
    checkAssembledQuantities(t, RRaw, RRawDense, totalDOFs);

    auto& KDense = denseFlatAssembler.getMatrix(req);
    auto& K      = sparseFlatAssembler.getMatrix(req);
    auto& RDense = denseFlatAssembler.getVector(req);
    auto& R      = sparseFlatAssembler.getVector(req);
    checkAssembledQuantities(t, K, KDense, totalDOFs);
    checkAssembledQuantities(t, R, RDense, totalDOFs);

    const Eigen::Index fixedDOFs = dirichletValues.fixedDOFsize();
    int boundaryNodes            = (elementsPerDirection[0] * Dune::power(2, ref) + 1) * 2 +
                        (elementsPerDirection[1] * Dune::power(2, ref) + 1) * 2 - 4;
    if constexpr (Ikarus::Concepts::LagrangeNodeOfOrder<
                      std::remove_cvref_t<decltype(fes[0].localView().tree().child(0))>, 2>)
      boundaryNodes *= 2;
    t.check(2 * boundaryNodes == fixedDOFs)
        << "Boundary DOFs (" << 2 * boundaryNodes << ") is not equal to Fixed DOFs (" << fixedDOFs << ")";

    /// check if full matrices and full vectors are correct after applying boundary conditions
    t.check(std::ranges::count(KDense.reshaped(), 1) == fixedDOFs) << "Correct number of ones in matrix";
    t.check(std::ranges::count(RDense, 0) == fixedDOFs) << "Correct number of zeros in vector";

    for (auto i = 0U; i < totalDOFs; ++i) {
      if (dirichletValues.isConstrained(i)) {
        for (auto j = 0U; j < totalDOFs; ++j) {
          t.check(Dune::FloatCmp::eq(KDense(i, j), static_cast<double>(i == j)))
              << std::string(i == j ? "" : "Off-") + "Diagonal term is " << KDense(i, j) << " at (" << i << "," << j
              << ")";
        }
        t.check(Dune::FloatCmp::eq(RDense(i), 0.0)) << "Vector component is " << RDense(i) << " at " << i;
      } else {
        for (auto j = 0U; j < totalDOFs; ++j)
          if (not(dirichletValues.isConstrained(j)))
            t.check(Dune::FloatCmp::eq(KDense(i, j), KRawDense(i, j)))
                << "KDense and KRawDense have the same entries corresponding to unconstrained DOFs";
        t.check(Dune::FloatCmp::eq(RDense(i), RRawDense(i)))
            << "RDense and RRawDense have the same entries corresponding to unconstrained DOFs";
      }
    }

    auto& KRedDense = denseFlatAssembler.getReducedMatrix(req);
    auto& KRed      = sparseFlatAssembler.getReducedMatrix(req);
    auto& RRedDense = denseFlatAssembler.getReducedVector(req);
    auto& RRed      = sparseFlatAssembler.getReducedVector(req);
    checkAssembledQuantities(t, KRed, KRedDense, totalDOFs - fixedDOFs);
    checkAssembledQuantities(t, RRed, RRedDense, totalDOFs - fixedDOFs);

    /// check if reduced matrices and reduced vectors are correct after applying boundary conditions
    Eigen::Index r = 0U;
    for (auto i = 0U; i < totalDOFs; ++i) {
      if (not(dirichletValues.isConstrained(i))) {
        Eigen::Index c = 0U;
        for (auto j = 0U; j < totalDOFs; ++j)
          if (not(dirichletValues.isConstrained(j))) {
            t.check(Dune::FloatCmp::eq(KRawDense(i, j), KRedDense(r, c)))
                << "Matrix components are " << KRawDense(i, j) << " and " << KRedDense(r, c) << " at (" << r << "," << c
                << ")";
            c++;
          }
        t.check(Dune::FloatCmp::eq(RRawDense(i), RRedDense(r)))
            << "Vector components are " << RRawDense(i) << " and " << RRedDense(r) << " at " << r;
        r++;
      }
    }

    grid->globalRefine(1);
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis  = power<2>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower2Basis = power<2>(lagrange<2>(), FlatInterleaved());

  t.subTest(SimpleAssemblersTest(firstOrderLagrangePrePower2Basis));
  t.subTest(SimpleAssemblersTest(secondOrderLagrangePrePower2Basis));
  return t.exit();
}
