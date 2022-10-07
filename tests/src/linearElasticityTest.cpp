#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
using Dune::TestSuite;

#include "common.hh"
#include "testHelpers.hh"

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/mechanics/linearElastic.hh>
#include <ikarus/finiteElements/mechanics/enhancedAssumedStrains.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/utils/algorithms.hh>

#include <vector>
#include <string>
template <int gridDim>
auto testLinearElasticityElement() {
  TestSuite t(std::string("testLinearElasticityElement"));

  using Grid = Dune::UGGrid<gridDim>;
  Dune::GridFactory<Grid> gridFactory;
  std::vector<Dune::FieldVector<double,gridDim>> corners;

  CornerFactory<gridDim>::construct(corners,4);

  for (auto& corner : corners) {
    gridFactory.insertVertex(corner);
  }

  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1,2,3});
  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();


  using namespace Ikarus;
//
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved()));
//
//  auto localView = basis.localView();

  auto volumeLoad = []([[maybe_unused]] auto& globalCoord, auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };
  const double youngsModulus = 1000;
  const double poissonsRatio = 1000;
  auto element = gridView.template begin<0>();
  Ikarus::EnhancedAssumedStrains<LinearElastic<decltype(basis)>> fe(basis, *element, youngsModulus, poissonsRatio);

    Eigen::VectorXd d;
    d.setZero(basis.size());
    double lambda = 0.0;

      Ikarus::FErequirements requirements = FErequirementsBuilder()
                                       .insertGlobalSolution(Ikarus::FESolutions::displacement, d)
                                       .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                                       .addAffordance(Ikarus::elastoStatics)
                                       .build();


Eigen::VectorXd forces;
Eigen::MatrixXd stiffnessmatrix;
      fe.calculateVector(requirements,forces);
      fe.calculateMatrix(requirements,stiffnessmatrix);


      std::cout<<"forces: "<<forces<<std::endl;
      std::cout<<"stiffnessmatrix: "<<stiffnessmatrix<<std::endl;
//
//  std::vector<bool> dirichletFlags(basis.size(), false);
//
//  Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
//    if (std::abs(intersection.geometry().center()[1]) < 1e-8) {
//      dirichletFlags[localView.index(localIndex)[0]] = true;
//    }
//  });
//
//  auto sparseAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);
//
//  Eigen::VectorXd d;
//  d.setZero(basis.size());
//  double lambda = 0.0;
//
//  auto residualFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
//    Ikarus::FErequirements req = FErequirementsBuilder()
//                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
//                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
//                                     .addAffordance(Ikarus::VectorAffordances::forces)
//                                     .build();
//    return sparseAssembler.getVector(req);
//  };
//
//  auto KFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
//    Ikarus::FErequirements req = FErequirementsBuilder()
//                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
//                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
//                                     .addAffordance(Ikarus::MatrixAffordances::stiffness)
//                                     .build();
//    return sparseAssembler.getMatrix(req);
//  };
//
//  auto energyFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
//    Ikarus::FErequirements req = FErequirementsBuilder()
//                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
//                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
//                                     .addAffordance(Ikarus::ScalarAffordances::mechanicalPotentialEnergy)
//                                     .build();
//    return sparseAssembler.getScalar(req);
//  };
//
//  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, KFunction),
//                                            parameter(d, lambda));
//
//  const double gradTol = 1e-8;
//
//  auto tr = Ikarus::makeTrustRegion(nonLinOp);
//  tr->setup({.verbosity = 1,
//             .maxiter   = 1000,
//             .grad_tol  = gradTol,
//             .corr_tol  = 1e-16,  // everything should converge to the gradient tolerance
//             .useRand   = false,
//             .rho_reg   = 1e8,
//             .Delta0    = 1});
//
//  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 2);
//  vtkWriter->setFileNamePrefix("Test2Dsolid");
//  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);
//
//  auto lc = Ikarus::LoadControl(tr, 1, {0, 2000});
//  lc.subscribeAll(vtkWriter);
//  const auto controlInfo = lc.run();
//  nonLinOp.template update<0>();
//  const auto maxDisp = std::ranges::max(d);
//  const double energyExpected
//      = (std::is_same_v<Grid, Grids::Yasp>)
//            ? -1.4809559783564966e+03
//            : ((std::is_same_v<Grid, Grids::Alu>) ? -1.4842107484533601e+03
//                                                  : /* std::is_same_v<Grid, Grids::Iga> */ -8.1142552237939071e+02);
//  const double maxDispExpected
//      = (std::is_same_v<Grid, Grids::Yasp>)
//            ? 0.786567027108437
//            : ((std::is_same_v<Grid, Grids::Alu>) ? 0.78426066482258983
//                                                  : /* std::is_same_v<Grid, Grids::Iga> */ 0.615624125459537153);
//
//  t.check(Dune::FloatCmp::eq(energyExpected, nonLinOp.value()), "energyExpected == nonLinOp.value()");
//  std::stringstream str;
//
//  t.check(std::abs(maxDispExpected - maxDisp) < 1e-12, "maxDispExpected-maxDisp");
//
//  nonLinOp.template update<1>();
//  t.check(controlInfo.success, "Successful result");
//  t.check(gradTol >= nonLinOp.derivative().norm(), "Gradient Tolerance should be larger than actual tolerance");
  return t;
}

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  t.subTest(testLinearElasticityElement<2>());
  t.subTest(testLinearElasticityElement<3>());

  return t.exit();
}