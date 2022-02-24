//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>

#include <dune/common/float_cmp.hh>
#include <dune/common/indices.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/onedgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/LocalBasis/localBasis.h"
#include "ikarus/basis/basishelper.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/FiniteElements/FEPolicies.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Solver/LinearSolver/LinearSolver.h>

template <typename Basis>
struct OneDimensionalTrussElement {
  using LocalView         = typename Basis::LocalView;
  using GlobalIndex       = typename LocalView::MultiIndex;
  using FERequirementType = Ikarus::FErequirements<Eigen::VectorXd>;
  OneDimensionalTrussElement(const Basis& basis, const typename LocalView::Element& element,
                             double p_EA) {
    // constructor
  }

  // Assembler requires calculateMatrix function
  [[nodiscard]] Eigen::Matrix2d calculateMatrix(const FERequirementType& requirements) const {
    if (requirements.matrixAffordances == Ikarus::stiffness)
      return stiffnessMatrix_;
    else
      throw std::runtime_error("The requested matrix type is not implemented;");
  }

  // Assembler requires globalIndices function
  [[nodiscard]] std::vector<GlobalIndex> globalIndices() const {
    return globalIndices_;
  }


private:
  Eigen::Matrix2d stiffnessMatrix_;
  std::vector<GlobalIndex> globalIndices_;
};

void exampleTrussElement(){
  constexpr int griddim = 1;
  const double L = 1;
  const double EA = 1.0;
  const int numElements = 10;
  const int numGP = 1;

  // create a grid and grid view
  Dune::OneDGrid grid(numElements, 0, L);
  auto gridView = grid.leafGridView();
  draw(gridView);

  // create a basis and a local view
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, lagrange<1>());
  auto localView = basis.localView();

  // initialize the global stiffness matrix
  auto numDofs = basis.size();
  Eigen::MatrixXd K_Glob = Eigen::MatrixXd::Zero(numDofs, numDofs);

  for (auto& ele : elements(gridView)){
    localView.bind(ele);
    auto& fe = localView.tree().finiteElement();
    Ikarus::LocalBasis localBasis(fe.localBasis());

    // define quadrature rule
    const auto& rule = Dune::QuadratureRules<double, 1>::
        rule(ele.type(),numGP,Dune::QuadratureType::GaussLegendre);

    // integration point loop
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(2,2);
    for (auto& gp : rule) {
      auto detJ = ele.geometry().integrationElement(0.0);
      Eigen::VectorXd dNdxi = Eigen::VectorXd::Zero(2);
      localBasis.evaluateJacobian(gp.position(),dNdxi);
      Eigen::VectorXd B = dNdxi/detJ;
      K += EA*B*B.transpose()*detJ*gp.weight();
    }

    // Adding local stiffness to the global stiffness
    for (auto i = 0U; i < localView.size(); ++i)
      for (auto j = 0U; j < localView.size(); ++j)
        K_Glob(localView.index(i)[0], localView.index(j)[0]) += K(i,j);
  }

  // external force (1 one the right end)
  Eigen::VectorXd F_ExtGlobRed = Eigen::VectorXd::Zero(numDofs-1);
  F_ExtGlobRed(numDofs-2) = 1.0;

  // reduce stiffness matrix (fix left end, i.e. fix first dof)
  const Eigen::MatrixXd K_GlobRed =
      K_Glob(Eigen::seq(1,Eigen::last),Eigen::seq(1,Eigen::last));

  // solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);
  linSolver.factorize(K_GlobRed);
  const Eigen::VectorXd D_GlobRed = linSolver.solve(F_ExtGlobRed);

  std::cout << "Displacement: " << D_GlobRed(Eigen::last) << "; exact solution: " << L/EA << "\n";

  }



void exampleTrussElementWithAssembler(){
  const double L = 1;
  const double EA = 1.0;
  const int numElements = 10;

  // grid
  Dune::OneDGrid grid(numElements, 0, L);
  auto gridView = grid.leafGridView();
  draw(gridView);

  // basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, lagrange<1>());
  auto numDofs = basis.size();

  // create vector of truss elements

  // fix system at the left end

  // create assembler and assemble stiffness matrix

  // external force (1 one the right end)
  Eigen::VectorXd F_ExtGlobRed = Eigen::VectorXd::Zero(numDofs);
  F_ExtGlobRed(numDofs-1) = 1.0;

  // solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);
  linSolver.factorize(K_GlobRed);
  const Eigen::VectorXd D_GlobRed = linSolver.solve(F_ExtGlobRed);

  std::cout << "Displacement: " << D_GlobRed(Eigen::last) << "; exact solution: " << L/EA << "\n";


}


  int main(){
  //exampleTrussElement();
  exampleTrussElementWithAssembler();
}