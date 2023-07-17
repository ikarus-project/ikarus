#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <dune/localfefunctions/eigenDuneTransformations.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
void getNodalBaseSytemMatrices(
    const auto&basisI,
    auto& directors
)
{
  auto basis=basisI.untouched();
  const auto gridView = basis.gridView();
  const int numCP= basis.size();
  Eigen::MatrixX3d A_ijGPN( numCP,3);
  A_ijGPN.setZero();
  Eigen::Matrix<double, 9, 1> A1A2A3;

  Eigen::MatrixXd NTN_ele;
  Eigen::MatrixX3d AijGPN_ele;
  auto localView= basis.localView();
  const auto _num_node = localView.size();
  AijGPN_ele.resize(_num_node, 3);
  AijGPN_ele.setZero();
  NTN_ele.resize(_num_node, _num_node);
  NTN_ele.setZero();

  Eigen::VectorXd N;
  std::vector<Eigen::Triplet<double>> tripletList;
  for (auto& ele: elements(gridView)) {
    NTN_ele.setZero();
    AijGPN_ele.setZero();
    localView.bind(ele);
    const int order=2 * localView.tree().finiteElement().localBasis().order();
    const auto localBasis = Dune::CachedLocalBasis(localView.tree().finiteElement().localBasis());
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(ele.type(), order);
    const auto geo= ele.geometry();
    for (auto gp: rule) {
      localBasis.evaluateFunction(gp.position(),N);

      const auto J = Dune::toEigen(geo.jacobianTransposed(gp.position()));
      const auto A3= (J.row(0).cross(J.row(1))).normalized().eval();

      for (int i = 0; i < _num_node; ++i)
        AijGPN_ele.row(i) += A3 * N[i];
    NTN_ele += N * N.transpose();
    }


    for (int i = 0; i < _num_node; ++i) ///Werte in Globales Gleichungs system einbauen
    {
      const auto I = localView.tree().localIndex(i);
      A_ijGPN.row(I) += AijGPN_ele.row(i);
      for (int j = 0; j < _num_node; ++j) {
        const auto J = localView.tree().localIndex(j);
        tripletList.push_back(Eigen::Triplet<double>(I,
                                                     J,
                                                     NTN_ele(i, j)));
      }
    }

  }
  Eigen::SimplicialLDLT<Eigen::SparseMatrix < double>, Eigen::Upper > solver;
  Eigen::SparseMatrix<double> NTN(numCP, numCP);
  NTN.setFromTriplets(tripletList.begin(), tripletList.end());
  solver.compute(NTN);
  Eigen::MatrixXd A_ijCP=solver.solve(A_ijGPN);
  std::cout<<A_ijCP<<std::endl;
  for (int i = 0; i < numCP; ++i) {
    directors[i].setValue(A_ijCP.row(i).transpose());
  }


}
