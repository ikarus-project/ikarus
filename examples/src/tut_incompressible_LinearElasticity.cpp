//
// Created by Alex on 21.07.2021.
//

#include <config.h>

#include <autodiff/forward/dual/dual.hpp>
#include <matplot/matplot.h>

#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/localFunctions/standardLocalFunction.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/utils/algorithms.hh>
#include <ikarus/utils/linearAlgebraTypedefs.hh>

using namespace Ikarus;
using namespace Dune::Indices;
template <typename Basis>
struct Solid : Ikarus::AutoDiffFEClean<Solid<Basis>, Basis> {
  using BaseAD = Ikarus::AutoDiffFEClean<Solid<Basis>, Basis>;
  friend BaseAD;
  using LocalView         = typename Basis::LocalView;
  using FERequirementType = typename BaseAD::FERequirementType;
  using Traits            = TraitsFromLocalView<LocalView>;
  Solid(const Basis& basis, const typename LocalView::Element& element, double emod, double nu)
      : BaseAD(basis, element), localView_{basis.localView()}, emod_{emod}, nu_{nu} {
    localView_.bind(element);

    mu_       = emod_ / (2 * (1 + nu_));
    lambdaMat = convertLameConstants({.emodul = emod_, .nu = nu_}).toLamesFirstParameter();
  }

  using GlobalIndex = typename LocalView::MultiIndex;
  void globalIndices(std::vector<GlobalIndex>& globalIndices) const {
    const auto& feDisp = localView_.tree().child(_0, 0).finiteElement();
    for (size_t i = 0; i < feDisp.size(); ++i) {
      for (int j = 0; j < Traits::worlddim; ++j) {
        globalIndices.push_back(localView_.index((localView_.tree().child(_0, j).localIndex(i))));
      }
    }
    const auto& fePressure = localView_.tree().child(_1).finiteElement();
    for (size_t i = 0; i < fePressure.size(); ++i) {
      globalIndices.push_back(localView_.index((localView_.tree().child(_1).localIndex(i))));
    }
  }

private:
  template <class ScalarType>
  [[nodiscard]] ScalarType calculateScalarImpl(const FERequirementType& par,
                                               const Eigen::VectorX<ScalarType>& dx) const {
    const auto& d      = par.getSolution(Ikarus::FESolutions::displacement);
    const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
    Eigen::VectorX<ScalarType> localDisp(localView_.size());
    localDisp.setZero();
    auto& displacementNode = localView_.tree().child(_0, 0);
    auto& pressureNode     = localView_.tree().child(_1);
    const auto& feDisp     = displacementNode.finiteElement();
    const auto& fePressure = pressureNode.finiteElement();
    Eigen::Matrix<ScalarType, Traits::dimension, Eigen::Dynamic> disp;
    disp.setZero(Eigen::NoChange, feDisp.size());
    for (auto i = 0U; i < feDisp.size(); ++i)
      for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
        disp.col(i)(k2) = dx[i * 2 + k2] + d[localView_.index(localView_.tree().child(_0, k2).localIndex(i))[0]];

    Eigen::Vector<ScalarType, Eigen::Dynamic> pN;
    pN.setZero(fePressure.size());
    for (auto i = 0U; i < fePressure.size(); ++i)
      pN[i] = dx[Traits::mydim * feDisp.size() + i] + d[localView_.index(localView_.tree().child(_1).localIndex(i))[0]];

    ScalarType energy = 0.0;

    const int order  = 2 * (feDisp.localBasis().order());
    const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
    const auto geo   = localView_.element().geometry();
    Ikarus::LocalBasis localBasisDisp(feDisp.localBasis());
    Ikarus::LocalBasis localBasisPressure(fePressure.localBasis());
    Eigen::Matrix<double, Eigen::Dynamic, Traits::mydim> dNdisp;
    Eigen::VectorXd Ndisp;
    Eigen::VectorXd Npressure;
    for (auto&& gp : rule) {
      const auto J = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
      localBasisDisp.evaluateFunctionAndJacobian(gp.position(), Ndisp, dNdisp);
      localBasisPressure.evaluateFunction(gp.position(), Npressure);
      const Eigen::Vector<double, Traits::worlddim> X = toEigenVector(geo.global(gp.position()));
      Eigen::Vector<ScalarType, Traits::worlddim> x   = X;
      for (auto i = 0U; i < feDisp.size(); ++i)
        x += disp.col(i) * Ndisp[i];

      ScalarType pressure = pN.dot(Npressure);

      const auto gradu      = (disp * dNdisp.template cast<ScalarType>() * J.inverse()).eval();
      const auto symgradu   = sym(gradu);
      const ScalarType divU = gradu.trace();

      Eigen::Vector<double, Traits::worlddim> fext;
      fext.setZero();
      fext[1] = lambda;
      fext[0] = 0 * lambda;

      energy += (0.5 * (2 * mu_ * symgradu.squaredNorm() - 1 / lambdaMat * Dune::power(pressure, 2)) + pressure * divU
                 - x.dot(fext))
                * geo.integrationElement(gp.position()) * gp.weight();  // plane strain for 2D
    }
    return energy;
  }

private:
  LocalView localView_;
  double emod_;
  double nu_;
  double mu_;
  double lambdaMat;
};

int main(int argc, char** argv) {
  using namespace Dune::Functions;
  /// Construct grid
  Dune::MPIHelper::instance(argc, argv);
  using namespace Ikarus;
  constexpr int gridDim = 2;
  //  /// ALUGrid Example
  //  using Grid = Dune::ALUGrid<gridDim, 2, Dune::cube, Dune::nonconforming>;
  //  auto grid  = Dune::GmshReader<Grid>::read("../../tests/src/testFiles/unstructuredQuadscoarse.msh", false);
  //  grid->globalRefine(1);
  //  auto gridView = grid->leafGridView();
  using Grid        = Dune::YaspGrid<gridDim>;
  const double L    = 1;
  const double h    = 1;
  const size_t elex = 20;
  const size_t eley = 20;

  Dune::FieldVector<double, 2> bbox = {L, h};
  std::array<int, 2> eles           = {elex, eley};
  auto grid                         = std::make_shared<Grid>(bbox, eles);
  auto gridView                     = grid->leafGridView();
  //  draw(gridView);

  using namespace Dune::Functions::BasisFactory;
  /// Construct basis
  auto basis
      = makeBasis(gridView, composite(power<2>(lagrange<1>(), FlatInterleaved()), lagrange<0>(), FlatLexicographic()));

  /// Create finite elements
  const double Emod = 2.1e1;
  const double nu   = 0.5;
  std::vector<Solid<decltype(basis)>> fes;
  for (auto& ele : elements(gridView))
    fes.emplace_back(basis, ele, Emod, nu);

  /// Collect dirichlet nodes
  std::vector<bool> dirichletFlags(basis.size(), false);
  forEachBoundaryDOF(subspaceBasis(basis, _0), [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = true;
  });

  /// Create assembler
  auto sparseFlatAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);
  auto denseFlatAssembler  = DenseFlatAssembler(basis, fes, dirichletFlags);

  /// Create non-linear operator
  double lambda = 0;
  Eigen::VectorXd d;
  d.setZero(basis.size());

  auto fintFunction = [&](auto&& lambdaLocal, auto&& dLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder().setSolution(Ikarus::FESolutions::displacement,dLocal).setParameter(Ikarus::FEParameter::loadfactor, lambdaLocal).setAffordance(Ikarus::VectorAffordances::forces).build();
    return denseFlatAssembler.getReducedVector(req);
  };
  auto KFunction = [&](auto&& lambdaLocal, auto&& dLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder().setSolution(Ikarus::FESolutions::displacement,dLocal).setParameter(Ikarus::FEParameter::loadfactor, lambdaLocal).setAffordance(Ikarus::MatrixAffordances::stiffness).build();
    return sparseFlatAssembler.getReducedMatrix(req);
  };

  auto K = KFunction(1.0, d);
  auto R = fintFunction(1.0, d);
  Eigen::SparseLU<decltype(K)> ld;
  ld.compute(K);
  if (ld.info() != Eigen::Success) assert(false && "Failed Compute");

  d -= denseFlatAssembler.createFullVector(ld.solve(R));
  if (ld.info() != Eigen::Success) assert(false && "Failed Solve");

  /// Postprocess
  auto disp
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(subspaceBasis(basis, _0), d);
  auto pressure = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _1), d);
  Dune::VTKWriter vtkWriter(gridView, Dune::VTK::nonconforming);
  vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.addVertexData(pressure, Dune::VTK::FieldInfo("pressure", Dune::VTK::FieldInfo::Type::scalar, 1));

  vtkWriter.write("TestIncompressible");
}