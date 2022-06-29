//
// Created by ac136645 on 6/15/2022.
//
#include <config.h>

#include <vector>

#include <dune/fufem/boundarypatch.hh>
#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/finiteElements/mechanics/Q1LinearElastic.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>

using namespace Ikarus;
using namespace Dune::Indices;

template <typename Basis>
class Q1LinearElasticAD : public DisplacementFE<Basis>,
                          public Ikarus::AutoDiffFE<Q1LinearElasticAD<Basis>, Basis> {
public:
  using BaseDisp = DisplacementFE<Basis>;  // Handles globalIndices function
  using BaseAD   = AutoDiffFE<Q1LinearElasticAD<Basis>, Basis>;
  using BaseAD::size;
  using GlobalIndex = typename DisplacementFE<Basis>::GlobalIndex;
  friend BaseAD;
  using FERequirementType = FErequirements<Eigen::VectorXd>;
  using LocalView         = typename Basis::LocalView;
  using GridView         = typename Basis::GridView;

  template <typename VolumeLoad, typename NeumannBoundaryLoad>
  Q1LinearElasticAD(Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu, const BoundaryPatch<GridView> * neumannBoundary,
                    const NeumannBoundaryLoad& neumannBoundaryLoad,
                    const VolumeLoad& p_volumeLoad)
      : BaseDisp(globalBasis, element),
        BaseAD(globalBasis, element),
        localView_{globalBasis.localView()},
        volumeLoad(p_volumeLoad),
        neumannBoundaryLoad_{neumannBoundaryLoad},
        neumannBoundary_{neumannBoundary},
        emod_{emod},
        nu_{nu} {
    localView_.bind(element);
    const int order = 2 * (localView_.tree().child(0).finiteElement().localBasis().order());
    localBasis      = Ikarus::LocalBasis(localView_.tree().child(0).finiteElement().localBasis());
    localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),
                    bindDerivatives(0, 1));
  }

  using Traits = TraitsFromLocalView<LocalView>;

private:
  template <class ScalarType>
  ScalarType calculateScalarImpl(const FERequirementType& req, Eigen::VectorX<ScalarType>& dx) const {
    const auto& d      = req.getSolution(Ikarus::FESolutions::displacement);
    const auto& lambda = req.getParameter(Ikarus::FEParameter::loadfactor);

    auto& first_child = localView_.tree().child(0);
    const auto& fe    = first_child.finiteElement();
    Dune::BlockVector<Ikarus::RealTuple<ScalarType, Traits::dimension>> disp(fe.size());

    for (auto i = 0U; i < fe.size(); ++i)
      for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
        disp[i][k2] = dx[i * 2 + k2] + d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];

    ScalarType energy = 0.0;
    const int order   = 2 * (fe.localBasis().order());
    const auto& rule  = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
    Eigen::Matrix3<ScalarType> C = planeStressLinearElasticMaterialTangent(emod_,nu_);
    const auto geo = localView_.element().geometry();
    Ikarus::StandardLocalFunction uFunction(localBasis, disp);
    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const auto Jinv = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().inverse().eval();
      const auto u    = uFunction.evaluateFunction(gpIndex);
      const auto H
          = uFunction.evaluateDerivative(gpIndex, wrt(DerivativeDirections::spatialAll), transformWith(Jinv));
      const auto E      = (0.5 * (H.transpose() + H)).eval();
      const auto EVoigt = toVoigt(E);

      Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigenVector(gp.position()), lambda);
      energy += (0.5 * EVoigt.dot(C * EVoigt) - u.dot(fext)) * geo.integrationElement(gp.position()) * gp.weight();
    }

    //line or surface loads, i.e. neumann boundary
    if (not neumannBoundary_) return energy;

    auto element = localView_.element();
    for (auto &&intersection : intersections(neumannBoundary_->gridView(), element))
    {
      if (not neumannBoundary_ or not neumannBoundary_->contains(intersection)) continue;

      const auto &quadLine = Dune::QuadratureRules<double, Traits::mydim-1>::rule(intersection.type(), order);

      for (const auto &curQuad : quadLine)
      {
        // Local position of the quadrature point
        const Dune::FieldVector<double, Traits::mydim> &quadPos = intersection.geometryInInside().global(curQuad.position());

        const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

        // The value of the local function
        const auto u    = uFunction.evaluateFunction(quadPos);

        // Value of the Neumann data at the current position
        auto neumannValue = neumannBoundaryLoad_(toEigenVector(intersection.geometry().global(curQuad.position())),lambda);

        energy -= neumannValue.dot(u) * curQuad.weight() * integrationElement;
      }
    }

    return energy;
  }

  LocalView localView_;
  Ikarus::LocalBasis<
      std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
      localBasis;
  //TODO: write as optional
  std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                        const double&)>
      volumeLoad;
  std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                        const double&)>
      neumannBoundaryLoad_;
  const BoundaryPatch<GridView> *neumannBoundary_;
  double emod_;
  double nu_;
};



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
    T0 << J11*J11 , J12*J12 , J11*J12 ,
          J21*J21 , J22*J22 , J21*J22 ,
      2.0*J11*J21 , 2.0*J12*J22 , J21*J12 + J11*J22;
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
  const auto Dinv = D.inverse();
  for (int i = 0; i < 4; ++i) {
    const auto& Li = L.block<4,2>(0,i*2);
    for (int j = 0; j < 4; ++j) {
      const auto& Lj = L.block<4,2>(0,j*2);
      K.block<2,2>(i*2,j*2)-= Li.transpose()*Dinv*Lj;
    }
  }
  return  K;
}


int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
  constexpr int gridDim = 2;
  const double E  = 1.0;
  const double nu = 1.0/3.0;

  double lambdaLoad = 1;

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

  // clamp left-hand side
  std::vector<bool> dirichletFlags(basis.size(),false);
  forEachBoundaryDOF(basis, [&](auto &&localIndex, auto &&localView, auto &&intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = true;
  });

  std::vector<Q1LinearElasticAD<decltype(basis)>> fesAD;
  std::vector<Ikarus::Q1LinearElastic<decltype(basis)>> fes;
  auto volumeLoad = [](auto& globalCoord, auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb*0;
    fext[0] = lamb*0;
    return fext;
  };

  auto neumannBoundaryLoad = [&](auto& globalCoord, auto& lamb) {
    Eigen::Vector2d F  = Eigen::Vector2d::Zero();
    F[1] = lamb/16.0;
    return F;
  };

  std::string lambdaNeumannVertices = std::string("lambda x: ( x[0]>47.9999 )");
  Python::start();
  Python::Reference main = Python::import("__main__");
  Python::run("import math");

  Python::runStream()
      << std::endl << "import sys"
      << std::endl << "import os"
      << std::endl;

  const auto& indexSet = gridView.indexSet();

  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  auto pythonNeumannVertices = Python::make_function<bool>(Python::evaluate(lambdaNeumannVertices));

  for (auto &&vertex: vertices(gridView))
  {
    bool isNeumann = pythonNeumannVertices(vertex.geometry().corner(0));
    neumannVertices[indexSet.index(vertex)] = isNeumann;
  }

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  for (auto& element : elements(gridView)) {
    auto localView = basis.localView();
    localView.bind(element);
    Q1E4Stiffness(localView,planeStressLinearElasticMaterialTangent(E,nu));
    fesAD.emplace_back(basis, element, E, nu, &neumannBoundary, neumannBoundaryLoad, volumeLoad);
    fes.emplace_back(basis, element, E, nu, &neumannBoundary, neumannBoundaryLoad, volumeLoad);
  }

  auto sparseAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);
  auto sparseAssemblerAD = SparseFlatAssembler(basis, fesAD, dirichletFlags);

  auto KFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                     .build();
    return sparseAssembler.getMatrix(req);
  };

  auto KFunctionAD = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                     .build();
    return sparseAssemblerAD.getMatrix(req);
  };


  auto residualFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::VectorAffordances::forces)
                                     .build();
    return sparseAssembler.getVector(req);
  };

  auto residualFunctionAD = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::VectorAffordances::forces)
                                     .build();
    return sparseAssemblerAD.getVector(req);
  };


  auto energyFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::ScalarAffordances::mechanicalPotentialEnergy)
                                     .build();
    return sparseAssembler.getScalar(req);
  };

  Eigen::VectorXd D_Glob = Eigen::VectorXd::Zero(basis.size());


  if(residualFunctionAD(D_Glob,lambdaLoad).isApprox(residualFunction(D_Glob,lambdaLoad)))
    std::cout<<"Coinciding external forces :)"<<std::endl;

  if(KFunctionAD(D_Glob,lambdaLoad).isApprox(KFunction(D_Glob,lambdaLoad)))
    std::cout<<"Coinciding stiffness :)"<<std::endl;


  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, KFunction),
                                                     parameter(D_Glob, lambdaLoad));


  const auto K = nonLinOp.secondDerivative();
  const auto Fext = -nonLinOp.derivative();

//  // solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_SimplicialLDLT);
  linSolver.compute(K);

  linSolver.solve(D_Glob, Fext);

  const auto deformedEnergy = nonLinOp.value();

  nonLinOp.update<0>();
  std::cout << "Energy after: " << nonLinOp.value() << std::endl;

  // Postprocess
  auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis, D_Glob);
  Dune::VTKWriter vtkWriter(gridView, Dune::VTK::conforming);
  vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.write("Cook_Membrane");
}
