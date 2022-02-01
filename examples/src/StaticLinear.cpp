#include <cmath>
#include <matplot/matplot.h>
#include <dune/alugrid/grid.hh>
#include <dune/grid/yaspgrid.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Assembler/SimpleAssemblers.h"
#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/FEManager/DefaultFEManager.h"
#include "ikarus/FiniteElements/ElasticityFE.h"
#include "ikarus/FiniteElements/NonLinearElasticityFE.h"
#include "ikarus/LinearAlgebra/DirichletConditionManager.h"
#include "ikarus/utils/Observer/controlLogger.h"
#include "ikarus/utils/Observer/gridDrawerObserver.h"
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include <ikarus/FiniteElements/ForceLoad.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include "ikarus/FiniteElements/NonLinearElasticityFEwithBasis.h"
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
template <typename Basis>
class DenseAssemblerFromBasis {
public:
  explicit DenseAssemblerFromBasis(const Basis& basis) : basis_{&basis} {}

  Eigen::MatrixXd& getMatrix(const Eigen::VectorXd& displacement, const double& lambda) {
    return getMatrixImpl(displacement, lambda);
  }

  Eigen::VectorXd& getVector(const Eigen::VectorXd& displacement, const double& lambda) {
    return getVectorImpl(displacement, lambda);
  }

private:
  Eigen::MatrixXd& getMatrixImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    mat.setZero(basis_->size(), basis_->size());
    auto localView = basis_->localView();
    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)> fe(localView, 1000, 0.0);
      auto matLoc = fe.calculateMatrix(displacement, lambda);
      auto first_child = localView.tree().child(0);
      for (auto i = 0U; i < localView.size(); ++i)
        for (auto j = 0U; j < localView.size(); ++j) {
          //          std::cout<<"LocalIndices: "<<localView.index(i)<<" "<<localView.index(j)<<" "<<i<<" "<<j<<std::endl;

          mat(localView.index(i)[0], localView.index(j)[0]) += matLoc(i, j);
        }
    }
    auto mat2 = mat;
    //    std::cout<<mat2<<std::endl;
    //    Eigen::FullPivLU<Eigen::MatrixXd> lu(mat2);
    //    std::cout<<"FullTRank: "<<lu.rank()<<std::endl;
    localView.unbind();
    for (auto i = 0U; i < 4; ++i)
      mat.col(i).setZero();
    for (auto i = 0U; i < 4; ++i)
      mat.row(i).setZero();
    for (auto i = 0U; i < 4; ++i)
      mat(i, i) = 1;
    //    std::cout << mat << std::endl;
    //    std::cout << "============================================" << std::endl;
    return mat;
  }

  Eigen::VectorXd& getVectorImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    vec.setZero(basis_->size());
    auto localView = basis_->localView();

    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      auto first_child = localView.tree().child(0);
      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)> fe(localView, 1000, 0.0);
      auto vecLocal = fe.calculateVector(displacement, lambda);
      for (auto i = 0U; i < localView.size(); ++i)
        vec(localView.index(i)[0]) += vecLocal(i);
    }
    for (auto i = 0U; i < 4; ++i)
      vec[i] = 0;
    return vec;
  }
  Basis const* basis_;
  Eigen::MatrixXd mat{};
  Eigen::VectorXd vec{};
};



int main() {
  constexpr int gridDim = 2;
//  using Grid            = Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>;
//  auto grid             = Dune::GmshReader<Grid>::read("../../tests/src/testFiles/unstructuredTrianglesfine.msh", false);
//  using GridView        = typename Grid::LeafGridView;

    using namespace Ikarus::Grid;
    using Grid = Dune::YaspGrid<gridDim>;
    const double L    = 1;
    const double h    = 1;
    const size_t elex = 10;
    const size_t eley = 10;

    Dune::FieldVector<double, 2> bbox = {L, h};
    std::array<int, 2> eles           = {elex, eley};
    auto grid                         = std::make_shared<Grid>(bbox, eles);

  using GridView    = typename Grid::LeafGridView;
  GridView gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved()));
  std::cout << "This gridview cotains: " << std::endl;
  std::cout << gridView.size(2) << " vertices" << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " elements" << std::endl;

  draw(gridView);
  std::cout << "1" << std::endl;

  auto denseAssembler = DenseAssemblerFromBasis(basis);

  Eigen::VectorXd d(basis.size());
  d.setZero();
  double lambda = 0.0;


  Dune::SubsamplingVTKWriter<GridView> vtkWriter(gridView, Dune::refinementLevels(2));
  //  Dune::VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(d, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, gridDim));
  vtkWriter.write("TestDuneBasis_0");

  double fac = 100;
  for (int ls = 1; ls < 10; ++ls) {
    Eigen::FullPivLU<Eigen::MatrixXd> lu;
    for (int i = 0; i < 20; ++i) {
      const auto& K = denseAssembler.getMatrix(d, lambda);
      lu.compute(K);
      const auto& r = denseAssembler.getVector(d, lambda);
      const auto dd = lu.solve(r);
      d-= dd;
      std::cout<<"Rnorm: "<<r.norm()<<"    "<<"dnorm: "<<dd.norm()<<"    "<<"Rank: "<<lu.rank()<<" Dofs: "<<lu.rows()<<std::endl;
      if(r.norm()<1e-8) break;
    }
    std::vector<Dune::FieldVector<double,2>> dv;
    dv.reserve(d.size()/2);
    for (auto pos=0U,i = 0U; i < dv.size(); ++i) {
      dv[i]= d(pos);
      dv[i+1]= d(pos+1);
      pos+=2;
    }
    std::cout<<d.transpose()<<std::endl;
    auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(basis,d);
    Dune::SubsamplingVTKWriter<GridView> vtkWriterI(gridView, Dune::refinementLevels(2));
    //    Dune::VTKWriter<GridView> vtkWriterI(gridView);
    vtkWriterI.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, gridDim));
    vtkWriterI.write("TestDuneBasis_" + std::to_string(ls));
    lambda+=fac;
  }
  return 0;
}