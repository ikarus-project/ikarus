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
    // create local view and bind it to an element
    LocalView localView = basis.localView();
    localView.bind(element);

    // calculate and store stiffness matrix
    double detJ = element.geometry().integrationElement(0.0);
    calculateStiffnessMatrix(localView, detJ, p_EA);

    // store glocal indices
    calculateGlobalIndices(localView);
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

  // function to calculate the stiffness matrix
  void calculateStiffnessMatrix(LocalView localView, double detJ, double EA) {
    auto& fe = localView.tree().finiteElement();
    Ikarus::LocalBasis localBasis(fe.localBasis());
    const auto& rule = Dune::QuadratureRules<double, 1>::
        rule(Dune::GeometryTypes::simplex(1),1,Dune::QuadratureType::GaussLegendre);

    Eigen::Vector2d dNdxi;
    stiffnessMatrix_.setZero();

    for (auto& gp : rule) {
      localBasis.evaluateJacobian(gp.position(),dNdxi);
      Eigen::Vector2d B = dNdxi/detJ;
      stiffnessMatrix_ += EA*B*B.transpose()*detJ*gp.weight();
    }
  }

  // function to get the global indices from the local view
  void calculateGlobalIndices(LocalView localView) {
    const auto& fe = localView.tree().finiteElement();
    for (size_t i = 0; i < fe.size(); ++i)
      globalIndices_.push_back(localView.index(localView.tree().localIndex(i)));
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

  // create a basis and a local view

  // initialize the global stiffness matrix

  // loop over all elements

    // bind basis to element and create local basis

    // define quadrature rule

    // integration point loop


    // Adding local stiffness to the global stiffness


  // external force (1 one the right end)

  // reduce stiffness matrix (fix left end, i.e. fix first dof)

  // solve the linear system

  // output the displacement

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
  exampleTrussElement();
  exampleTrussElementWithAssembler();
}