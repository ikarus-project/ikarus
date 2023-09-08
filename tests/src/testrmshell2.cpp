// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>
#include <ranges>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/nurbsgrid.hh>
#include <ikarus/finiteElements/feBases/autodiffFE.hh>

#include "spdlog/spdlog.h"
#include "testCommon.hh"

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/mechanics/fesettings.hh>
#include <ikarus/finiteElements/mechanics/nonlinearrmshell.hh>
#include <ikarus/finiteElements/mechanics/rmlstressBased.hh>
#include <ikarus/io/resultFunction.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>

#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template<typename CASStrain>
auto testMembraneStrain(const auto& localView,const auto& d)
{
  TestSuite t;
  const auto element = localView.element();
   auto geoR = element.geometry();

   const auto cps= geoR.impl().controlPoints().directGetAll();
  for (auto cp:cps) {
    std::cout<<cp<<std::endl;
  }
  auto geo = std::make_shared<const decltype(geoR)>(   element.geometry());
  auto &first_child = localView.tree().child(0);
  const auto &fe = first_child.finiteElement();
  const auto localBasis = Dune::CachedLocalBasis(fe.localBasis());
  Dune::BlockVector<Dune::RealTuple<double, 3>> disp(fe.size());
  for (auto i = 0U; i < disp.size(); ++i)
    for (auto k2 = 0U; k2 < 3; ++k2)
      disp[i][k2] = d[localView.index(localView.tree().child(k2).localIndex(i))[0]];

  Dune::FieldVector<double, 2> gpPos={0.34,0.57};
  auto f = [&](auto&& dx){
    Dune::BlockVector<Dune::RealTuple<autodiff::dual, 3>> dispD(fe.size());
      for (auto i = 0U; i < dispD.size(); ++i)
        for (auto k2 = 0U; k2 < 3; ++k2)
          dispD[i][k2]
              = dx[i*3 + k2] + disp[i][k2];

    Dune::StandardLocalFunction uFunction(localBasis, dispD, geo);
    CASStrain strain;
    return strain.value(gpPos,geoR,uFunction);
  };

  Dune::StandardLocalFunction uFunction(localBasis, disp, geo);

  Eigen::VectorXdual dx(localView.size());
  Eigen::VectorXdual g(localView.size());
  dx.setZero();
  Eigen::MatrixXd h=jacobian(f, autodiff::wrt(dx), at(dx), g);

  CASStrain strain;
  Eigen::MatrixXd hR(3,localView.size());
  for (int i = 0; i < fe.size(); ++i) {
    hR.block<3,3>(0,3*i)=strain.derivative(gpPos,Eigen::Matrix<double, 2, 3>(),double(),geoR,uFunction,localBasis,i);
  }

  const double tol = 1e-10;
  std::cout<<"h\n"<<h<<std::endl;
  std::cout<<"hR\n"<<hR<<std::endl;

  t.check(h.isApprox(hR, tol))<<
                              "Mismatch between the  matrices obtained from explicit implementation and the one based on "
                              "automatic differentiation: h"<<"\n" << h <<"\n hR \n"<< hR<<"\n h-hR \n"<< h-hR;

  return t;
}
//
//
template<template<typename> typename ShellElement>
auto checkFEByAutoDiff(std::string filename) {
  TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation of Kirchhoff-Love shell "+filename);

  constexpr auto dimworld        = 3;
  const std::array<int, 2> order = {1,1};

//  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0,  1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

//  const std::vector<std::vector<ControlPoint>> controlPoints
//      = {
//          {{.p = {0, 0.0, 0}, .w = 1}, {.p = {0, 1, 0}, .w = 1}, {.p = {0, 2, 0}, .w = 1}},
//          {{.p = {5, 0.0, 0}, .w = 1}, {.p = {5, 1, 0}, .w = 1}, {.p = {5, 2, 0}, .w = 1}},
//          {{.p = {10, 0.0, 0}, .w = 1}, {.p = {10, 1, 0}, .w = 1}, {.p = {10, 2, 0}, .w = 1}}};
//
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {0, 1, 0}, .w = 1}}, {{.p = {12, 0, 0}, .w = 1}, {.p = {12, 1, 0}, .w = 1}}};

  std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

  Dune::IGA::NURBSPatchData<2, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
//  for (int i = 0; i < 2; ++i)
//    patchData = degreeElevate(patchData, i, 1);
  auto grid = std::make_shared<Grid>(patchData);
  grid->globalDegreeElevate(1);
  grid->globalRefine(0);
  grid->globalDegreeElevate(1);


  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto scalarMidSurfBasis = nurbs();
//  auto scalaarDirectorBasis = Imp::NurbsPreBasisFactory<2, 3>();
  auto scalaarDirectorBasis = Imp::NurbsPreBasisFactory<2, 3>(gridView.impl().lowerOrderPatchData(1));
  auto basis       = Ikarus::makeBasis(gridView, composite(power<3>(scalarMidSurfBasis),power<2>(scalaarDirectorBasis)));
  auto element     = gridView.template begin<0>();
  auto nDOF        = basis.flat().size();
  auto localView        = basis.flat().localView();

  localView.bind(*element);
  auto nDOFPerEle        = localView.size();
//  Eigen::VectorXd dT;
//  dT.setZero(nDOF);
//  t.subTest(testMembraneStrain<Ikarus::CASMembraneStrain<Ikarus::CASAnsatzFunction>>(localView,dT));
//  std::cout<<"========================="<<std::endl;
//  t.subTest(testMembraneStrain<Ikarus::CASMembraneStrain<Ikarus::CASAnsatzFunctionANS>>(localView,dT));
  const double tol = 1e-10;

  auto volumeLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    std::array<Eigen::Vector<double,3>,2> vLoad;
    auto& fext = vLoad[0];
    auto& mext = vLoad[1];
    fext.setZero();
//    fext[0] = lamb;
//    fext[1] = 0.01 * lamb;
    fext[1] = 2*lamb;
    mext.setZero();
    return vLoad;
  };

  auto neumannBoundaryLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    std::array<Eigen::Vector<double,3>,2> vLoad;
    auto& fext = vLoad[0];
    auto& mext = vLoad[1];
    fext.setZero();
//    fext[0] = lamb;
//    fext[1] = 0.01 * lamb;
    fext[2] = 4*lamb;
    mext.setZero();
//    mext[0]=lamb;
    return vLoad;
  };

  /// We artificially apply a Neumann load on the complete boundary
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);
  std::array<std::string,3> directorFunctions;
  directorFunctions[0]="NFE";
  directorFunctions[1]="PBFE";
  directorFunctions[2]="GFE";
  for (int i = 0; i < 1; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 1; ++k) {
        Ikarus::FESettings feSettings;
        feSettings.addOrAssign("youngs_modulus", 1000.0);
        feSettings.addOrAssign("poissons_ratio", 0.0);
        feSettings.addOrAssign("thickness", 0.1);
        feSettings.addOrAssign("membraneStrainFlag", i);
        feSettings.addOrAssign("transverseShearStrainFlag", k);
        feSettings.addOrAssign("directorFunction", directorFunctions[j]);
        feSettings.addOrAssign("secondOrderBending", true);
        using Basis = decltype(basis);
        //  KLSHELL fe(basis, *element, feSettings);
        using MidSurfaceVector   = Dune::BlockVector<Dune::RealTuple<double, 3>>;
        using DirectorVector     = Dune::BlockVector<Dune::UnitVector<double, 3>>;
        using MultiTypeVector    = Dune::TupleVector<MidSurfaceVector, DirectorVector>;
        using MultiTypeVectorRaw = Dune::TupleVector<Dune::BlockVector<Dune::FieldVector<double, 3>>,
                                                     Dune::BlockVector<Dune::FieldVector<double, 3>>>;
        using namespace Dune::Functions::BasisFactory;
        auto blockedmidSurfaceBasis = Dune::Functions::subspaceBasis(basis.untouched(), Dune::Indices::_0);

        MidSurfaceVector mBlocked(basis.untouched().size({Dune::Indices::_0}));
        //    auto refCoords = [](auto v){ return Dune::FieldVector<double,3>();};
        //    Functions::interpolate(blockedmidSurfaceBasis, mBlocked, refCoords);
        //    auto deformationPowerBasis = makeBasis(gridView,power<3>(nurbs()));

        for (auto& msingle : mBlocked) {
          //      msingle.setValue(Eigen::Vector<double, 3>::Zero());
          msingle.setValue(Eigen::Vector<double, 3>::Random());
        }

        DirectorVector dBlocked(basis.untouched().size({Dune::Indices::_1}));
        for (auto& dsingle : dBlocked) {
          dsingle.setValue(Eigen::Vector<double, 3>::UnitZ() + 0.1 * Eigen::Vector<double, 3>::Random());
          //      dsingle.setValue(Eigen::Vector<double, 3>::UnitZ());
        }

        const MultiTypeVector x0(mBlocked, dBlocked);
        MultiTypeVector x = x0;
        //  ShellElement<Basis> fe(basis, *element, feSettings,x0, volumeLoad, &neumannBoundary, neumannBoundaryLoad);
        std::vector<ShellElement<Basis>> fes;
        for (auto& ele : elements(gridView))
          fes.emplace_back(basis, ele, feSettings, x0, volumeLoad, &neumannBoundary, neumannBoundaryLoad);
        using AutoDiffBasedFE = Ikarus::AutoDiffFE<ShellElement<Basis>, true>;
        //  AutoDiffBasedFE feAutoDiff(fe);

        //  Eigen::VectorXd d;
        //  d.setRandom(nDOF);
        // d.setZero(nDOF);
        //    auto basis3D       = Ikarus::makeBasis(gridView, composite(power<3>(scalarMidSurfBasis, BlockedInterleaved()),power<3>(scalaarDirectorBasis, BlockedInterleaved()),
        //                                                               BlockedLexicographic{}));
        //    auto blockedmidSurfaceBasis2 = Dune::Functions::subspaceBasis(basis3D.untouched(),Dune::Indices::_0);
        //    auto blockeddirectorBasis2 = Dune::Functions::subspaceBasis(basis3D.untouched(),Dune::Indices::_1);
        //    auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(blockedmidSurfaceBasis2,  x); auto director = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(blockeddirectorBasis2,  x);
        //  Dune::SubsamplingVTKWriter vtkWriter(gridView,Dune::refinementLevels(0));
        //
        //    vtkWriter.addVertexData(disp, {"displacements", Dune::VTK::FieldInfo::Type::vector, 3});
        //    vtkWriter.addVertexData(director, {"director", Dune::VTK::FieldInfo::Type::vector, 3});
        //  vtkWriter.write(filename+ std::to_string(i));

        Ikarus::DirichletValues dirichletValues(basis.flat());
        auto sparseAssembler = Ikarus::SparseFlatAssembler(fes, dirichletValues);

        //  auto localDisp=localFunction(disp);
        //  localDisp.bind(*element);

        double lambda = 7.3;

        auto req = Ikarus::FErequirements<std::reference_wrapper<MultiTypeVector>>()
                       .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                       .addAffordance(Ikarus::AffordanceCollections::elastoStatics);

        auto fvLambda = [&](auto&& d_) -> auto {
          req.insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, d_);
          return sparseAssembler.getScalar(req);
        };

        auto dfvLambda = [&](auto&& d_) -> auto& {
          req.insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, d_);
          return sparseAssembler.getVector(req);
        };
        auto ddfvLambda = [&](auto&& d_) -> auto& {
          req.insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, d_);
          return sparseAssembler.getMatrix(req);
        };
        auto nonLinOp
            = Ikarus::NonLinearOperator(Ikarus::functions(fvLambda, dfvLambda, ddfvLambda), Ikarus::parameter(x));

        nonLinOp.updateAll();

        auto uF = std::function([&](MultiTypeVector& multiTypeVector, const Eigen::VectorXd& d_) {
          //      auto dFull = sparseAssembler.createFullVector(d_);
          //    std::cout<<"dFull"<<std::endl;
          //    std::cout<<dFull<<std::endl;
          using Ikarus::operator+=;
          multiTypeVector += d_;
        });

        t.check(Ikarus::checkGradient(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}, uF))
            << "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
            << "automatic differentiation with membraneStrainFlag: " << i << " and director func "
            << directorFunctions[j] <<" TransversShearStrain: "<<k<< "\n  The residual is \n " << dfvLambda(x);
        t.check(Ikarus::checkHessian(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}, uF))
            << "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
            << "automatic differentiation with membraneStrainFlag: " << i << " and director func "
            << directorFunctions[j] <<" TransversShearStrain: "<<k<< "\n The Hessian is\n"
            << ddfvLambda(x);
        //    auto subOp = nonLinOp.template subOperator<1, 2>();
        //    t.check(Ikarus::checkJacobian(subOp, {.draw = false, .writeSlopeStatementIfFailed = true},uF))<<"Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "<<
        //        "automatic differentiation of the gradient with membraneStrainFlag: "<<i<<" and director func "<<directorFunctions[j];

        //  Eigen::MatrixXd K, KAutoDiff;
        //  K.setZero(nDOFPerEle, nDOFPerEle);
        //  KAutoDiff.setZero(nDOFPerEle, nDOFPerEle);
        //
        //  Eigen::VectorXd R, RAutoDiff;
        //  R.setZero(nDOFPerEle);
        //  RAutoDiff.setZero(nDOFPerEle);
        //
        //  fe.calculateVector(req, R);
        //  feAutoDiff.calculateVector(req, RAutoDiff);
        //    fe.calculateMatrix(req, K);
        //    feAutoDiff.calculateMatrix(req, KAutoDiff);
        //
        //  t.check(K.isApprox(KAutoDiff, tol),"K Check"+filename)<<
        //      "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
        //      "automatic differentiation with membraneStrainFlag: "<<i<<" and director func "<<directorFunctions[j]<<"\n" << K <<"\n KAutoDiff \n"<< KAutoDiff<<"\n K-KAutoDiff \n"<< K-KAutoDiff;
        //
        //  t.check(R.isApprox(RAutoDiff, tol),"R Check"+filename)<<
        //      "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
        //      "automatic differentiation with membraneStrainFlag: "<<i<<" and director func "<<directorFunctions[j]<<"\n" << R <<"\n RAutoDiff \n"<< RAutoDiff<<"\n R-RAutoDiff \n"<< R-RAutoDiff;
        //
        //  t.check(Dune::FloatCmp::eq(fe.calculateScalar(req), feAutoDiff.calculateScalar(req), tol),"E Check"+filename)<<
        //    "Mismatch between the energies obtained from explicit implementation and the one based on "
        //    "automatic differentiation"<<"with membraneStrainFlag: "<<i<<" and director func "<<directorFunctions[j];
      }
    }
  }
  return t;
}




template <typename Basis>
using RMSHELL = Ikarus::NonLinearRMshell<Basis>;
template <typename Basis>
using RMSHELLSB = Ikarus::StressBasedShellRM<RMSHELL<Basis>>;
#include <cfenv>
int main(int argc, char** argv) {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  Ikarus::init(argc, argv);
  //  const double E             = materialParameters.get<double>("E");
  //  const double nu            = materialParameters.get<double>("nu");

//  checkFEByAutoDiff<RMSHELL>("RMSHELL");
  std::cout<<"RMSHELLSB"<<std::endl;
  checkFEByAutoDiff<RMSHELLSB>("RMSHELLSB");

//  checkFEByAutoDiff<KLSHELLSB>("KLSHELLSB");
}
