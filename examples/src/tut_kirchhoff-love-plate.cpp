//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>

#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/nurbsgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/LocalBasis/localBasis.h"
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/FiniteElements/AutodiffFE.h>
#include <ikarus/utils/concepts.h>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include "ikarus/utils/Observer/controlVTKWriter.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/Controlroutines/LoadControl.h"

struct KirchhoffPlate {
  static constexpr double Emodul = 2.1e8;
  static constexpr double nu = 0.3;
  static constexpr double thickness = 0.1;

  static  Eigen::Matrix<double,3,3> constitutiveMatrix(double Emod, double nu, double thickness)
  {
    const double factor = Emod*Dune::power(thickness,3)/(12*(1-nu*nu));
    Eigen::Matrix<double,3,3> D;
    D.setZero();
    D(0,0) = 1;
    D(0,1) = D(1,0) = nu;
    D(1,1)  = 1;
    D(2,2)  = (1-nu)/2.0;
    return D;
  }

  template <typename LocalView, class Scalar>
  static Scalar calculateScalarImpl(const LocalView& localView, const Eigen::VectorXd& x,
                                    const Eigen::VectorX<Scalar>& dx, const double lambda) {
    const auto D = constitutiveMatrix(Emodul,nu,thickness);
    Scalar energy = 0.0;
    auto& ele     = localView.element();
    const auto X1 = Ikarus::toEigenVector(ele.geometry().corner(0));
    const auto X2 = Ikarus::toEigenVector(ele.geometry().corner(1));
    auto& fe = localView.tree().finiteElement();
    Eigen::VectorX<Scalar> wNodal;
    wNodal.setZero(fe.size());
    for (auto i = 0U; i < fe.size(); ++i)
      wNodal(i) = dx[localView.tree().localIndex(i)] + x[localView.index(localView.tree().localIndex(i))[0]];

    const auto& localBasis = fe.localBasis();

    const auto& rule = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * localBasis.order());
    for (auto& gp : rule) {
      std::vector<Dune::FieldMatrix<double, 1, 2>> dN_xi_eta;
      std::vector<Dune::FieldVector<double, 1>> dN_xixi;
      std::vector<Dune::FieldVector<double, 1>> dN_xieta;
      std::vector<Dune::FieldVector<double, 1>> dN_etaeta;
      std::vector<Dune::FieldVector<double, 1>> N_dune;
      Eigen::VectorXd N(fe.size());

      localBasis.evaluateJacobian(gp.position(), dN_xi_eta);
      localBasis.evaluateFunction(gp.position(), N_dune);
      std::ranges::copy(N_dune, N.begin());
      localBasis.partial({2, 0}, gp.position(), dN_xixi);
      localBasis.partial({1, 1}, gp.position(), dN_xieta);
      localBasis.partial({0, 2}, gp.position(), dN_etaeta);
      const auto Jinv = Ikarus::toEigenMatrix(ele.geometry().jacobianInverseTransposed(gp.position())).transpose().eval();
      Eigen::VectorXd dN_xx(fe.size());
      Eigen::VectorXd dN_yy(fe.size());
      Eigen::VectorXd dN_xy(fe.size());
      for (auto i = 0U; i < fe.size(); ++i) {
        dN_xx[i] = dN_xixi[i]*Jinv(0, 0) + dN_etaeta[i]*Jinv(0, 1);
        dN_yy[i] = dN_xixi[i]*Jinv(1, 1) + dN_etaeta[i]*Jinv(1, 0);
        dN_xy[i] = dN_xi_eta[i][0][0]*Jinv(0, 0)*Jinv(1, 0) + dN_xi_eta[i][0][1]*Jinv(1, 1)*Jinv(0, 1);
      }
      Eigen::Vector<Scalar,3> kappa;
      kappa(0) = dN_xx.dot(wNodal);
      kappa(1)=dN_yy.dot(wNodal);
      kappa(2) =dN_xy.dot(wNodal);
      Scalar w = N.dot(wNodal);

      energy+= (0.5 * kappa.dot(D*kappa)-w*lambda)  * ele.geometry().integrationElement(gp.position())* gp.weight();
    }


    return energy;
  }


};


template <typename Basis>
class DenseFlatAssembler {
 public:
  explicit DenseFlatAssembler(const Basis& basis, const std::vector<bool>& dirichFlags)
      : basis_{&basis}, dirichletFlags{&dirichFlags} {}

  Eigen::MatrixXd& getMatrix(const Eigen::VectorXd& displacement, const double& lambda) {
    return getMatrixImpl(displacement, lambda);
  }

  Eigen::VectorXd& getVector(const Eigen::VectorXd& displacement, const double& lambda) {
    return getVectorImpl(displacement, lambda);
  }

  double getScalar(const Eigen::VectorXd& displacement, const double& lambda) {
    return getScalarImpl(displacement, lambda);
  }

 private:
  Eigen::MatrixXd& getMatrixImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    mat.setZero(basis_->size(), basis_->size());
    auto localView = basis_->localView();
    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      auto matLoc      = Ikarus::AutoDiffFE<KirchhoffPlate>::calculateMatrix(localView, displacement,lambda);
      for (auto i = 0U; i < localView.size(); ++i)
        for (auto j = 0U; j < localView.size(); ++j) {
          mat(localView.index(i)[0], localView.index(j)[0]) += matLoc(i, j);
        }
    }
    for (auto i = 0U; i < basis_->size(); ++i)
      if (dirichletFlags->at(i)) mat.col(i).setZero();
    for (auto i = 0U; i < basis_->size(); ++i)
      if (dirichletFlags->at(i)) mat.row(i).setZero();
    for (auto i = 0U; i < basis_->size(); ++i)
      if (dirichletFlags->at(i)) mat(i, i) = 1;
    return mat;
  }

  Eigen::VectorXd& getVectorImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    vec.setZero(basis_->size());
    auto localView = basis_->localView();

    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      auto vecLocal = Ikarus::AutoDiffFE<KirchhoffPlate>::calculateVector(localView, displacement,lambda);
      for (auto i = 0U; i < localView.size(); ++i)
        vec(localView.index(i)[0]) += vecLocal(i);
    }
    for (auto i = 0U; i < basis_->size(); ++i) {
      if (dirichletFlags->at(i)) vec[i] = 0;
    }

    return vec;
  }

  double getScalarImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    double scalar = 0.0;
    vec.setZero(basis_->size());
    auto localView = basis_->localView();

    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      for (auto i = 0U; i < localView.size(); ++i)
        scalar += Ikarus::AutoDiffFE<KirchhoffPlate>::calculateScalar(localView, displacement,lambda);
    }

    return scalar;
  }
  Basis const* basis_;
  std::vector<bool> const* dirichletFlags;
  Eigen::MatrixXd mat{};
  Eigen::VectorXd vec{};
};





int main() {
  using namespace Ikarus;
  constexpr int griddim                                    = 2;
  constexpr int dimworld                                   = 2;
  const std::array<std::vector<double>, griddim> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointType;

  const double Lx = 1;
  const double Ly = 1;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {Lx, 0}, .w = 1}}, {{.p = {0, Ly}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, griddim> dimsize = {2, 2};

  auto controlNet = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<griddim, dimworld>;

  Dune::IGA::NURBSPatchData<griddim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = {1, 1};
  patchData.controlPoints = controlNet;
  patchData               = Dune::IGA::degreeElevate(patchData, 0, 1);
  patchData               = Dune::IGA::degreeElevate(patchData, 1, 1);
  Grid grid(patchData);
  grid.globalRefine(2);
  auto gridView = grid.leafGridView();
  draw(gridView);
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, gridView.getPreBasis());
  std::vector<double> w(basis.size());
  auto localView = basis.localView();

  auto dirichletPredicate = [](auto p) { return std::sin(p[0]); };
  std::vector<double> uhat(basis.size());

  std::vector<bool> dirichletFlags(basis.size());
  std::fill(dirichletFlags.begin(), dirichletFlags.end(), false);

  Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& index) { dirichletFlags[index]=true; });
  auto denseAssembler = DenseFlatAssembler(basis, dirichletFlags);

  Eigen::VectorXd d;
  d.setZero(basis.size());
  double lambda = 0.0;

  auto fintFunction   = [&](auto&& lambda, auto&& disp) -> auto& { return denseAssembler.getVector(disp, lambda); };
  auto KFunction      = [&](auto&& lambda, auto&& disp) -> auto& { return denseAssembler.getMatrix(disp, lambda); };
  auto energyFunction = [&](auto&& lambda, auto&& disp) -> auto { return denseAssembler.getScalar(disp, lambda); };

  auto nonLinOp  = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, fintFunction, KFunction),
                                             parameter(lambda, d));
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);

  auto nr                      = Ikarus::NewtonRaphson(nonLinOp.subOperator<1, 2>(), std::move(linSolver));
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 2);
  vtkWriter->setFileNamePrefix("TestKLplate");
  vtkWriter->setVertexSolutionName("displacement");
  nr.subscribeAll(nonLinearSolverObserver);

  auto lc = Ikarus::LoadControl(std::move(nr), 20, {0, 2000});

  lc.subscribeAll(vtkWriter);
  std::cout << "Energy before: " << nonLinOp.value() << std::endl;
  lc.run();
  nonLinOp.update<0>();
  std::cout << "Energy after: " << nonLinOp.value() << std::endl;



}