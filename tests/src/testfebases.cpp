// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/nedelecbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/grid/yaspgrid.hh>

using Dune::TestSuite;

#include <ikarus/finiteelements/febases.hh>
#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

/**
 * \brief LeafNodeVisitor helps to perform specific operations
 * at a leaf node, see dune/typetree/visitor.hh for more details
 *
 * \tparam LeafOp Type of the functor to perform operations at the leaf node
 */
template <class LeafOp = Dune::TypeTree::NoOp>
struct LeafNodeVisitor : Dune::TypeTree::TreeVisitor, Dune::TypeTree::DynamicTraversal
{
  explicit LeafNodeVisitor(const LeafOp& leafOp = {})
      : leafOp_(leafOp) {}

  template <class Node, class TreePath>
  void leaf(Node&& node, TreePath tp) const {
    leafOp_(node, tp);
  }

private:
  LeafOp leafOp_;
};

/**
 * \brief A helper function to get a FEBase.
 *
 * \tparam GridView Type of the grid view
 * \tparam PreBasis Type of the pre basis
 * *
 * \param gridView The DUNE grid view.
 * \param preBasis The pre basis to make a basis.
 *
 * \return FEBase corresponding to the gridView and preBasis.
 */
template <typename GridView, typename PreBasis>
auto getBasis(const GridView& gridView, const PreBasis& preBasis) {
  const auto basis      = Ikarus::makeBasis(gridView, preBasis);
  const auto& localView = basis.flat().localView();
  auto element          = elements(gridView).begin();
  return Ikarus::FEBase<decltype(basis)>(basis, *element);
}

/**
 * \brief A helper function to get global indices using LeafNodeVisitor.
 *
 * \tparam Basis Type of the FEBase
 *
 * \param basis The FE base to access localView functionalities.
 *
 * \return Output vector to store global indices.
 */
template <typename Basis>
auto leafNodeIndices(const Basis& basis) {
  std::vector<typename std::remove_cvref_t<decltype(basis.localView())>::MultiIndex> dofs;

  auto leafOpFunc = [&](auto&& node, auto&&) {
    const auto& fe = node.finiteElement();
    for (size_t i = 0; i < fe.size(); ++i) {
      dofs.push_back(basis.localView().index(node.localIndex(i)));
    }
  };

  auto visitor = LeafNodeVisitor(leafOpFunc);

  Dune::TypeTree::applyToTree(basis.localView().tree(), visitor);

  return dofs;
}

/**
 * \brief A test suite to test a set of bases.
 * \tparam numberOfBases Number of bases to be tested.
 * \tparam Bases Type of different bases.
 * \param expectedNumberOfChildren Expected number of children for each basis.
 * \param expectedSize Expected size for each basis.
 * \param bases A tuple of bases.
 * \param basisName Name of the types of bases.
 * \return Test suite object.
 */
template <size_t numberOfBases, typename... Bases>
auto BasisTest(const std::array<int, numberOfBases>& expectedNumberOfChildren,
               const std::array<int, numberOfBases>& expectedSize, const std::tuple<Bases...>& bases,
               const std::string& basisName) {
  static_assert(std::tuple_size_v<std::tuple<Bases...>> == numberOfBases, "Input size mismatch in BasisTest.");
  TestSuite t("BasisTest");
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<numberOfBases>()), [&](const auto i) {
    const auto basis          = std::get<i>(bases);
    const std::string message = " for i = " + std::to_string(i) + " -> " + basisName;
    using MultiIndex          = typename std::remove_cvref_t<decltype(basis.localView())>::MultiIndex;
    static_assert(Dune::IsIndexable<MultiIndex>(), "MultiIndex must support operator[]");
    checkScalars(t, basis.numberOfChildren(), expectedNumberOfChildren[i],
                 " Number of children is incorrect" + message);
    checkScalars(t, static_cast<int>(basis.size()), expectedSize[i], " Size is incorrect" + message);
    std::vector<MultiIndex> dofs;
    Ikarus::FEHelper::globalFlatIndices(basis.localView(), dofs);
    checkScalars(t, dofs.size(), basis.size(), " Size of dofs and basis is not equal" + message);
    t.check(!dofs.empty(), "dofs is empty" + message);
    auto leafNodeDofs = leafNodeIndices(basis);
    std::ranges::sort(dofs);
    std::ranges::sort(leafNodeDofs);
    t.check(dofs == leafNodeDofs, "dofs is not equal to leafNodeDofs after sorting" + message);
    t.check(dofs[0] == 0, "Smallest index contains a non-zero entry" + message);
    const bool hasDuplicates = std::adjacent_find(dofs.begin(), dofs.end()) == dofs.end();
    t.check(hasDuplicates) << "The sorted dofs vector has duplicates" + message;
  });
  return t;
}

auto FEBaseTest() {
  TestSuite t("FEBaseTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {1, 1};
  const auto grid                         = std::make_shared<Grid>(bbox, elementsPerDirection);

  const auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;

  const auto firstOrderLagrangePreBasis  = lagrange<1>();
  const auto secondOrderLagrangePreBasis = lagrange<2>();

  const auto firstOrderLagrangePowerPreBasis  = power<2>(firstOrderLagrangePreBasis);
  const auto secondOrderLagrangePowerPreBasis = power<5>(secondOrderLagrangePreBasis);

  const auto scalarScalarCompositePreBasis = composite(firstOrderLagrangePreBasis, secondOrderLagrangePreBasis);
  const auto scalarPowerCompositePreBasis  = composite(firstOrderLagrangePreBasis, secondOrderLagrangePowerPreBasis);

  const auto powerScalarCompositePreBasis = composite(firstOrderLagrangePowerPreBasis, secondOrderLagrangePreBasis);
  const auto powerPowerCompositePreBasis = composite(firstOrderLagrangePowerPreBasis, secondOrderLagrangePowerPreBasis);

  const auto scalarCompositePreBasis = composite(firstOrderLagrangePreBasis);
  const auto powerCompositePreBasis  = composite(firstOrderLagrangePowerPreBasis);

  const auto combinedPreBasis =
      composite(composite(scalarCompositePreBasis, powerPowerCompositePreBasis, firstOrderLagrangePreBasis),
                secondOrderLagrangePowerPreBasis);

  const auto compositePowerCombinedBasis = composite(power<10>(combinedPreBasis), secondOrderLagrangePreBasis);

  const auto lagrangeDGPowerPreBasis     = power<6>(lagrangeDG<3>());
  const auto nedelecScalarPreBasis       = nedelec<1, 1>();
  const auto raviartThomasScalarPreBasis = raviartThomas<2>();
  const auto specialCompositePreBasis =
      composite(lagrangeDGPowerPreBasis, nedelecScalarPreBasis, raviartThomasScalarPreBasis);

  /// Types of scalar bases
  const auto scalar1 = getBasis(gridView, firstOrderLagrangePreBasis);
  const auto scalar2 = getBasis(gridView, secondOrderLagrangePreBasis);

  /// Types of power bases
  const auto power1 = getBasis(gridView, firstOrderLagrangePowerPreBasis);
  const auto power2 = getBasis(gridView, secondOrderLagrangePowerPreBasis);

  /// Types of composite bases
  const auto composite1 = getBasis(gridView, scalarScalarCompositePreBasis);
  const auto composite2 = getBasis(gridView, scalarPowerCompositePreBasis);
  const auto composite3 = getBasis(gridView, powerScalarCompositePreBasis);
  const auto composite4 = getBasis(gridView, powerPowerCompositePreBasis);
  const auto composite5 = getBasis(gridView, scalarCompositePreBasis);
  const auto composite6 = getBasis(gridView, powerCompositePreBasis);
  const auto composite7 = getBasis(gridView, combinedPreBasis);
  const auto composite8 = getBasis(gridView, compositePowerCombinedBasis);

  /// Types of special bases
  const auto special1 = getBasis(gridView, lagrangeDGPowerPreBasis);
  const auto special2 = getBasis(gridView, nedelecScalarPreBasis);
  const auto special3 = getBasis(gridView, raviartThomasScalarPreBasis);
  const auto special4 = getBasis(gridView, specialCompositePreBasis);

  t.subTest(BasisTest<2>({0, 0}, {4, 9}, std::make_tuple(scalar1, scalar2), "Scalar basis"));
  t.subTest(BasisTest<2>({2, 5}, {8, 45}, std::make_tuple(power1, power2), "Power Basis"));
  t.subTest(BasisTest<6>({2, 2, 2, 2, 1, 1}, {13, 49, 17, 53, 4, 8},
                         std::make_tuple(composite1, composite2, composite3, composite4, composite5, composite6),
                         "Composite Basis"));
  t.subTest(BasisTest<1>({2}, {106}, std::make_tuple(composite7), "Combined Composite Basis"));
  t.subTest(BasisTest<1>({2}, {1069}, std::make_tuple(composite8), "Combined Composite and Power Basis"));
  t.subTest(BasisTest<1>({6}, {96}, std::make_tuple(special1), "Lagrange DG Basis"));
  t.subTest(BasisTest<1>({0}, {4}, std::make_tuple(special2), "Nedelec Basis"));
  t.subTest(BasisTest<1>({0}, {24}, std::make_tuple(special3), "Raviart Thomas Basis"));
  t.subTest(BasisTest<1>({3}, {124}, std::make_tuple(special4), "Special Composite Basis"));

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(FEBaseTest());
  return t.exit();
}
