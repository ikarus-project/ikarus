//
// Created by ac136645 on 6/15/2022.
//
#include <config.h>
#include <vector>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>


using namespace Ikarus;
using namespace Dune::Indices;

Eigen::MatrixXd Q1Stiffness(auto localView, const Eigen::Matrix3d& C) {

  using namespace Dune::Indices;
  constexpr int gridDim   = 2;
  auto element            = localView.element();
  auto geometry           = element.geometry();
  Eigen::MatrixXd K       = Eigen::MatrixXd::Zero(localView.size(), localView.size());

  const auto& localFiniteElement = localView.tree().child(_0).finiteElement();
  // Define the integration rule
  int order = 2;
  const auto& quadRule
      = Dune::QuadratureRules<double, gridDim>::rule(element.type(), order, Dune::QuadratureType::GaussLegendre);

  for (const auto& quadPoint : quadRule)
  {
    const auto quadPos = quadPoint.position();
    const auto jacobianinvT = geometry.jacobianInverseTransposed(quadPos); //J^{-1}.Transpose() in Dune = J^{-1}
    const auto detJ = geometry.integrationElement(quadPos); //determinant(J)

    std::vector<Dune::FieldMatrix<double,1,gridDim>> referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos,referenceGradients);
    std::vector<Dune::FieldVector<double,gridDim>> gradients(referenceGradients.size());

    for (size_t i=0; i<gradients.size(); i++)
      jacobianinvT.mv(referenceGradients[i][0],gradients[i]);

    // setup B-operator
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, localView.size());
    Eigen::VectorXd dNdx = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd dNdy = Eigen::VectorXd::Zero(4);
    for(size_t i=0;i<4;i++)
    {
      dNdx[i] = gradients[i][0];
      dNdy[i] = gradients[i][1];
    }

    for(size_t i=0;i<4;i++)
    {
      B(0,i) = dNdx[i];
      B(2,i) = dNdy[i];

      B(1,i+4) = dNdy[i];
      B(2,i+4) = dNdx[i];
    }
    // integration of stiffness matrix
    K += B.transpose() * C * B * detJ * quadPoint.weight();
  }

//  Eigen::EigenSolver<Eigen::MatrixXd> esk(K);
//  std::cout << "The eigenvalues of K are:" << std::endl << esk.eigenvalues() << std::endl;

  return K;
}

Eigen::MatrixXd Q1E4Stiffness(auto localView, const Eigen::Matrix3d& C) {

  using namespace Dune::Indices;
  constexpr int gridDim   = 2;
  auto element            = localView.element();
  auto geometry           = element.geometry();
  Eigen::MatrixXd K       = Eigen::MatrixXd::Zero(localView.size(), localView.size());
  Eigen::MatrixXd L       = Eigen::MatrixXd::Zero(4, localView.size());
  Eigen::MatrixXd D       = Eigen::MatrixXd::Zero(4, 4);

  const auto& localFiniteElement = localView.tree().child(_0).finiteElement();
  // Define the integration rule
  int order = 2;
  const auto& quadRule
      = Dune::QuadratureRules<double, gridDim>::rule(element.type(), order, Dune::QuadratureType::GaussLegendre);

  for (const auto& quadPoint : quadRule)
  {
    const auto quadPos = quadPoint.position();
    const auto jacobianinvT = geometry.jacobianInverseTransposed(quadPos); //J^{-1}.Transpose() in Dune = J^{-1}
    const auto detJ = geometry.integrationElement(quadPos); //determinant(J)

    Dune::FieldVector<double,2> quadPos0;
    quadPos0[0] = 0.5; // Center of the Element in Domain [0,1]
    quadPos0[1] = 0.5; // Center of the Element in Domain [0,1]

    const auto jacobianinvT0 = geometry.jacobianInverseTransposed(quadPos0); //J^{-1}.Transpose() in Dune = J^{-1}
    const auto detJ0 = geometry.integrationElement(quadPos0); //determinant(J)

    Eigen::MatrixXd jaco_it = Eigen::MatrixXd::Zero(2,2);
    for (size_t i=0;i<2;i++)
      for (size_t j=0;j<2;j++)
        jaco_it(i,j)=jacobianinvT0[i][j];

    auto jaco = (jaco_it).inverse();
    auto J11 = jaco(0,0);
    auto J12 = jaco(0,1);
    auto J21 = jaco(1,0);
    auto J22 = jaco(1,1);

    Eigen::Matrix3d T0;
    T0 << J11*J11 , J12*J12 , J11*J12 , J21*J21 , J22*J22 , J21*J22 , 2.0*J11*J21 , 2.0*J12*J22 , J21*J12 + J11*J22;
    auto T0_inv = T0.inverse() * (detJ0/detJ);

    std::vector<Dune::FieldMatrix<double,1,gridDim>> referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos,referenceGradients);
    std::vector<Dune::FieldVector<double,gridDim>> gradients(referenceGradients.size());

    for (size_t i=0; i<gradients.size(); i++)
      jacobianinvT.mv(referenceGradients[i][0],gradients[i]);

    // setup B-operator
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, localView.size());
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(3, 4);
    Eigen::VectorXd dNdx = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd dNdy = Eigen::VectorXd::Zero(4);
    for(size_t i=0;i<4;i++)
    {
      dNdx[i] = gradients[i][0];
      dNdy[i] = gradients[i][1];
    }

    for(size_t i=0;i<4;i++)
    {
      B(0,i) = dNdx[i];
      B(2,i) = dNdy[i];

      B(1,i+4) = dNdy[i];
      B(2,i+4) = dNdx[i];
    }

    // -(1/2) is added in the M-Matrix in order to fulfill the Orthogonality condition for Q1E4 element in domain [0,1]
    M(0,0) = quadPos[0]-0.5;
    M(1,1) = quadPos[1]-0.5;
    M(2,2) = quadPos[0]-0.5;
    M(2,3) = quadPos[1]-0.5;

    M = T0_inv * M ;

    // integration of stiffness matrix
    K += B.transpose() * C * B * detJ * quadPoint.weight();
    L += M.transpose() * C * B * detJ * quadPoint.weight();
    D += M.transpose() * C * M * detJ * quadPoint.weight();
  }

//  Eigen::EigenSolver<Eigen::MatrixXd> esk(K);
//  std::cout << "The eigenvalues of K are:" << std::endl << esk.eigenvalues() << std::endl;
//
//  Eigen::EigenSolver<Eigen::MatrixXd> eskeas(K - (L.transpose() * D.inverse() * L));
//  std::cout << "The eigenvalues of Keas are:" << std::endl << eskeas.eigenvalues() << std::endl;

  return K - (L.transpose() * D.inverse() * L);
}

Eigen::VectorXd get_line_load(auto localView, auto gridView, Eigen::Vector2d load)
{
  using namespace Dune::Indices;
  constexpr int gridDim   = 2;
  auto element            = localView.element();
  auto geometry           = element.geometry();
  Eigen::VectorXd fe       = Eigen::VectorXd::Zero(localView.size());
  int order = 2;

  const auto& localFiniteElement = localView.tree().child(_0).finiteElement();

  for (const auto& intersection : intersections(gridView, element))
  {
    if (intersection.geometry().center()[0]>47.99)
    {
      const auto &quadLine = Dune::QuadratureRules<double, gridDim-1>::rule(intersection.type(), order);
      for (const auto &curQuad : quadLine)
      {
        // Local position of the quadrature point
        const Dune::FieldVector<double, gridDim> &quadPos = intersection.geometryInInside().global(curQuad.position());
        const double detJ = intersection.geometry().integrationElement(curQuad.position()); //integrationElement evaluated on boundary
        std::vector<Dune::FieldVector<double,1>> shapeFunctionValues;
        localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);
        Eigen::Matrix2Xd N = Eigen::Matrix2Xd::Zero(2,8);
        for(size_t nn=0;nn<4;nn++)
        {
          N(0,nn) = shapeFunctionValues[nn];
          N(1,nn+4) = shapeFunctionValues[nn];
        }
        fe += N.transpose() * load * detJ * curQuad.weight();
      }
    }
  }
  return fe;
}



int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
  constexpr int gridDim = 2;
  const double E  = 1.0;
  const double nu = 1.0/3.0;
  const double En = E/(1.0-(nu*nu));
  const double Es = (1.0-nu)/2.0;
  Eigen::Vector2d F  = Eigen::Vector2d::Zero();
  F[1] = 1.0/16.0;

  Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
  C << En, En*nu, 0, En*nu, En, 0, 0,0,En*Es;

  using Grid = Dune::UGGrid<gridDim>;
  auto grid  = Dune::GmshReader<Grid>::read("../../tests/src/testFiles/cook.msh", false);
  grid->globalRefine(2);
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved()));

  // FlatInterleaved() -> <ux1,ux2,ux3,ux4,uy1,uy2,uy3,uy4>

  std::cout << "This gridview contains: " << std::endl;
  std::cout << gridView.size(2) << " vertices" << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " elements" << std::endl;
  std::cout << basis.size() << " Dofs" << std::endl;

  draw(gridView);
  auto localView = basis.localView();

  // global stiffness matrix and force vector
  auto numDofs              = basis.size();
  Eigen::VectorXd F_ExtGlob = Eigen::VectorXd::Zero(numDofs);
  Eigen::MatrixXd K_Glob    = Eigen::MatrixXd::Zero(numDofs, numDofs);

  for (auto& ele : elements(gridView)) {
    localView.bind(ele);

    // get local stiffness matrix and local load vector
    auto K_local = Q1E4Stiffness(localView, C);
    auto f_local = get_line_load(localView,gridView,F);
    // Adding local stiffness the global stiffness
    for (auto i = 0U; i < localView.size(); ++i)
    {
      for (auto j = 0U; j < localView.size(); ++j)
        K_Glob(localView.index(i)[0], localView.index(j)[0]) += K_local(i, j);
      F_ExtGlob(localView.index(i)[0]) += f_local(i);
    }

  }

  // clamp left-hand side
  Eigen::VectorXi dirichletFlags = Eigen::VectorXi::Zero(basis.size());
  forEachBoundaryDOF(subspaceBasis(basis, _0), [&](auto &&localIndex, auto &&localView, auto &&intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = 1;
  });
  forEachBoundaryDOF(subspaceBasis(basis, _1), [&](auto &&localIndex, auto &&localView, auto &&intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = 1;
  });

  for (size_t i=0;i<basis.size();i++)
  {
    if (dirichletFlags[i]==1)
    {
      K_Glob.col(i).setZero();
      K_Glob.row(i).setZero();
      K_Glob(i, i) = 1.0;
    }
  }

  // solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);
  linSolver.factorize(K_Glob);
  Eigen::VectorXd D_Glob = Eigen::VectorXd::Zero(numDofs);
  linSolver.solve(D_Glob, F_ExtGlob);

  // Postprocess
  auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis, D_Glob);
  Dune::VTKWriter vtkWriter(gridView, Dune::VTK::conforming);
  vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.write("Cook_Membrane");
}
