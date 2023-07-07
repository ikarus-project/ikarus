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
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/nurbsgrid.hh>
#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

#include "spdlog/spdlog.h"
#include "ikarus/io/vtkFunctionExtensions.hh"
#include "ikarus/io/composedGgridfuncMod.hh"

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/mechanics/fesettings.hh>
#include <ikarus/finiteElements/mechanics/nonlinearrmshell.hh>
#include <ikarus/finiteElements/mechanics/rmlstressBased.hh>
#include <ikarus/io/resultFunction.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/solver/nonLinearSolver/trustRegion.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>
#include <ikarus/io/shell3DDataCollector.hh>

using Dune::TestSuite;
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

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


template<template<typename> typename ShellElement>
auto checkFEByAutoDiff(std::string filename) {
  TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation of Kirchhoff-Love shell");

  constexpr auto dimworld        = 3;
  const std::array<int, 2> order = {2,2};

  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
//  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0,  1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {
          {{.p = {0, 0.0, 0}, .w = 1}, {.p = {0, 1, 0}, .w = 1}, {.p = {0, 2, 0}, .w = 1}},
          {{.p = {5, 0.0, 0}, .w = 1}, {.p = {5, 1, 0}, .w = 1}, {.p = {5, 2, 0}, .w = 1}},
          {{.p = {10, 0.0, 0}, .w = 1}, {.p = {10, 1, 0}, .w = 1}, {.p = {10, 2, 0}, .w = 1}}};

//  const std::vector<std::vector<ControlPoint>> controlPoints
//      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {0, 1, 0}, .w = 1}}, {{.p = {12, 0, 0}, .w = 1}, {.p = {12, 1, 0}, .w = 1}}};

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

//  grid->globalRefine(2);
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto scalarMidSurfBasis = nurbs();
  auto scalaarDirectorBasis = nurbs();
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
  for (int i = 0; i < 3; ++i) {


  Ikarus::FESettings feSettings;
  feSettings.addOrAssign("youngs_modulus", 1000.0);
  feSettings.addOrAssign("poissons_ratio", 0.0);
  feSettings.addOrAssign("thickness", 0.1);
  feSettings.addOrAssign("simulationFlag", i);
  using Basis = decltype(basis);
//  KLSHELL fe(basis, *element, feSettings);
    using MidSurfaceVector = Dune::BlockVector<Dune::RealTuple<double, 3>>;
    using DirectorVector  = Dune::BlockVector<Dune::UnitVector<double, 3>>;
    using MultiTypeVector = Dune::TupleVector<MidSurfaceVector, DirectorVector>;
    using MultiTypeVectorRaw = Dune::TupleVector< Dune::BlockVector<Dune::FieldVector<double, 3>>, Dune::BlockVector<Dune::FieldVector<double, 3>>>;
    using namespace Dune::Functions::BasisFactory;
    auto blockedmidSurfaceBasis = Dune::Functions::subspaceBasis(basis.untouched(),Dune::Indices::_0);

    MidSurfaceVector mBlocked(basis.untouched().size({Dune::Indices::_0}));
//    auto refCoords = [](auto v){ return Dune::FieldVector<double,3>();};
//    Functions::interpolate(blockedmidSurfaceBasis, mBlocked, refCoords);
//    auto deformationPowerBasis = makeBasis(gridView,power<3>(nurbs()));


    for (auto &msingle : mBlocked) {
//      msingle.setValue(Eigen::Vector<double, 3>::Zero());
      msingle.setValue(Eigen::Vector<double, 3>::Random());
    }


    DirectorVector dBlocked(basis.untouched().size({Dune::Indices::_1}));
    for (auto &dsingle : dBlocked) {
      dsingle.setValue(Eigen::Vector<double, 3>::UnitZ()+0.1*Eigen::Vector<double, 3>::Random());
//      dsingle.setValue(Eigen::Vector<double, 3>::UnitZ());
    }

    const MultiTypeVector x0(mBlocked, dBlocked);
    MultiTypeVector x =x0;
  ShellElement<Basis> fe(basis, *element, feSettings,x0, volumeLoad, &neumannBoundary, neumannBoundaryLoad);
  using AutoDiffBasedFE = Ikarus::AutoDiffFE<ShellElement<Basis>, true>;
  AutoDiffBasedFE feAutoDiff(fe);

//  Eigen::VectorXd d;
//  d.setRandom(nDOF);
  //d.setZero(nDOF);
    auto basis3D       = Ikarus::makeBasis(gridView, composite(power<3>(scalarMidSurfBasis, BlockedInterleaved()),power<3>(scalaarDirectorBasis, BlockedInterleaved()),
                                                               BlockedLexicographic{}));
    auto blockedmidSurfaceBasis2 = Dune::Functions::subspaceBasis(basis3D.untouched(),Dune::Indices::_0);
    auto blockeddirectorBasis2 = Dune::Functions::subspaceBasis(basis3D.untouched(),Dune::Indices::_1);
    auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(blockedmidSurfaceBasis2,  x);
    auto director = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(blockeddirectorBasis2,  x);
  Dune::SubsamplingVTKWriter vtkWriter(gridView,Dune::refinementLevels(0));

    vtkWriter.addVertexData(disp, {"displacements", Dune::VTK::FieldInfo::Type::vector, 3});
    vtkWriter.addVertexData(director, {"director", Dune::VTK::FieldInfo::Type::vector, 3});
  vtkWriter.write(filename+ std::to_string(i));

//  auto localDisp=localFunction(disp);
//  localDisp.bind(*element);

  double lambda = 7.3;

    auto req = Ikarus::FErequirements<MultiTypeVector>()
        .insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, x)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
        .addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  Eigen::MatrixXd K, KAutoDiff;
  K.setZero(nDOFPerEle, nDOFPerEle);
  KAutoDiff.setZero(nDOFPerEle, nDOFPerEle);

  Eigen::VectorXd R, RAutoDiff;
  R.setZero(nDOFPerEle);
  RAutoDiff.setZero(nDOFPerEle);

  fe.calculateVector(req, R);
  feAutoDiff.calculateVector(req, RAutoDiff);
    fe.calculateMatrix(req, K);
    feAutoDiff.calculateMatrix(req, KAutoDiff);

  t.check(K.isApprox(KAutoDiff, tol),"K Check"+filename)<<
      "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
      "automatic differentiation with simulationFlag: "<<i<<"\n" << K <<"\n KAutoDiff \n"<< KAutoDiff<<"\n K-KAutoDiff \n"<< K-KAutoDiff;

  t.check(R.isApprox(RAutoDiff, tol),"R Check"+filename)<<
      "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
      "automatic differentiation with simulationFlag: "<<i<<"\n" << R <<"\n RAutoDiff \n"<< RAutoDiff<<"\n R-RAutoDiff \n"<< R-RAutoDiff;

  t.check(Dune::FloatCmp::eq(fe.calculateScalar(req), feAutoDiff.calculateScalar(req), tol),"E Check"+filename)<<
    "Mismatch between the energies obtained from explicit implementation and the one based on "
    "automatic differentiation"<<"with simulationFlag: "<<i;
  }
  return t;
}

auto NonLinearElasticityLoadControlNRandTRforRMShell() {
  TestSuite t("NonLinearElasticityLoadControlNRandTRforKLShell ");
  constexpr auto dimworld        = 3;
  const std::array<int, 2> order = {1, 1};

  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;
  const double b = 1;
  const double L = 12;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {0, b, 0}, .w = 1}}, {{.p = {L, 0, 0}, .w = 1}, {.p = {L, b, 0}, .w = 1}}};

  std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

  Dune::IGA::NURBSPatchData<2, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;


  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree("/tmp/Ikarus/tests/src/shell.parset",
                                         parameterSet);

  const auto E              = parameterSet.get<double>("E");
  const auto nu             = parameterSet.get<double>("nu");
  const auto thickness      = parameterSet.get<double>("thickness");
  const auto loadFactor     = parameterSet.get<double>("loadFactor");
  const auto simulationFlag = parameterSet.get<int>("simulationFlag");
  const auto refine         = parameterSet.get<int>("refine");
  const auto orderElevate         = parameterSet.get<std::array<int,2>>("orderElevate");
  const auto plotInPlaneRefine         = parameterSet.get<int>("plotInPlaneRefine");
  const auto loadSteps         = parameterSet.get<int>("loadSteps");

  auto grid = std::make_shared<Grid>(patchData);
  for (int i = 0; i < 2; ++i)
    grid->degreeElevateInDirection(i,orderElevate[i]);
  grid->globalRefineInDirection(0,refine);

  auto gridView = grid->leafGridView();



  using GridView = decltype(gridView);
  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;

  // const double E         = 100000;
  // const double nu        = 0.0;
  // const double thickness = 0.001;
  Ikarus::FESettings feSettings;
  feSettings.addOrAssign("youngs_modulus", E);
  feSettings.addOrAssign("poissons_ratio", nu);
  feSettings.addOrAssign("thickness", thickness);
  feSettings.addOrAssign("simulationFlag", simulationFlag);

  auto scalarMidSurfBasis = nurbs();
  auto scalaarDirectorBasis = nurbs();
  auto basis       = Ikarus::makeBasis(gridView, composite(power<3>(scalarMidSurfBasis, BlockedInterleaved()),power<2>(scalaarDirectorBasis, BlockedInterleaved()),
                                                           BlockedLexicographic{}));

  std::cout<<"Number of Elements: "<<gridView.size(0)<<std::endl;
  std::cout<<"Dofs: "<<basis.flat().dimension()<<std::endl;
//  auto volumeLoad = [thickness, loadFactor]([[maybe_unused]] auto& globalCoord, auto& lamb) {
//    Eigen::Vector3d fext;
//    fext.setZero();
//    //    fext[1]= 2 * Dune::power(thickness, 3) * lamb / 10;
//    fext[2] = 2 * Dune::power(thickness, 3) * lamb * loadFactor;
//    return fext;
//  };
  auto volumeLoad = [thickness, loadFactor]<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    std::array<Eigen::Vector<double,3>,2> vLoad;
    auto& fext = vLoad[0];
    auto& mext = vLoad[1];
    fext.setZero();
//    fext[0] = lamb;
//    fext[1] = 0.01 * lamb;
    fext[2] = 2 * Dune::power(thickness, 3) * lamb * loadFactor*0;
    mext.setZero();

    return vLoad;
  };

  auto boundaryLoad = [thickness, loadFactor,E,L]<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    std::array<Eigen::Vector<double,3>,2> vLoad;
    auto& fext = vLoad[0];
    auto& mext = vLoad[1];
    fext.setZero();
    //    fext[0] = lamb;
    //    fext[1] = 0.01 * lamb;
    fext[2] = 2 * Dune::power(thickness, 3) * lamb * loadFactor*0;
    mext.setZero();
    const double pi         = std::numbers::pi;
    mext[1]=2.0*pi*E*Dune::power(thickness, 3)/L/12.0*(1.0-thickness*thickness*pi*pi/(L*L))*lamb*loadFactor;
//    std::cout<<"vLoad[0]"<<std::endl;
//    std::cout<<vLoad[0]<<std::endl;
//    std::cout<<"vLoad[1]"<<std::endl;
//    std::cout<<vLoad[1]<<std::endl;
//    std::cout<<"loadFactor"<<std::endl;
//    std::cout<<loadFactor<<std::endl;
//    std::cout<<"lamb"<<std::endl;
//    std::cout<<lamb<<std::endl;
    return vLoad;
  };

  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  const auto& indexSet= gridView.indexSet();
  for(auto& vertex: vertices(gridView))
    if(vertex.geometry().center()[0]>12-1e-8)
      neumannVertices[indexSet.index(vertex)]=true;
  std::cout<<neumannVertices<<std::endl;
  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);


  using MidSurfaceVector = Dune::BlockVector<Dune::RealTuple<double, 3>>;
  using DirectorVector  = Dune::BlockVector<Dune::UnitVector<double, 3>>;
  using MultiTypeVector = Dune::TupleVector<MidSurfaceVector, DirectorVector>;
  using MultiTypeVectorRaw = Dune::TupleVector< Dune::BlockVector<Dune::FieldVector<double, 3>>, Dune::BlockVector<Dune::FieldVector<double, 3>>>;
  using namespace Dune::Functions::BasisFactory;
  auto blockedmidSurfaceBasis = Dune::Functions::subspaceBasis(basis.untouched(),Dune::Indices::_0);

  MidSurfaceVector mBlocked(basis.untouched().size({Dune::Indices::_0}));

  for (auto &msingle : mBlocked) {
    msingle.setValue(Eigen::Vector<double, 3>::Zero());
  }


  DirectorVector dBlocked(basis.untouched().size({Dune::Indices::_1}));
  for (auto &dsingle : dBlocked) {
    dsingle.setValue(Eigen::Vector<double, 3>::UnitZ());
  }

  const MultiTypeVector x0(mBlocked, dBlocked);
  MultiTypeVector x =x0;

  using ElementTypePRim = Ikarus::NonLinearRMshell<decltype(basis)>;
  using ElementTypeRaw = Ikarus::StressBasedShellRM<ElementTypePRim>;
//  using ElementType = Ikarus::AutoDiffFE<ElementTypePRim, true>;
  using ElementType = ElementTypeRaw;
  std::vector<ElementType> fes;

  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, feSettings,x0, volumeLoad,&neumannBoundary,boundaryLoad);

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());

  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
  });

//  dirichletValues.fixDOFs([&](auto& basis, auto&& dirichletFlags) {
//    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basis,Dune::Indices::_0),
//                                        [&](auto&& localIndex, auto&& localView, auto&& intersection) {
//                                          if (std::abs(intersection.geometry().center()[0]) < 1e-8)
//                                            dirichletFlags[localView.index(localIndex)] = true;
//                                        });
//  });
//
//  dirichletValues.fixDOFs([&](auto& basis, auto&& dirichletFlags) {
//    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basis,Dune::Indices::_0, 2),
//                                        [&](auto&& localIndex, auto&& localView, auto&& intersection) {
//                                          if (std::abs(intersection.geometry().center()[0]) > 10 - 1e-8)
//                                            dirichletFlags[localView.index(localIndex)] = true;
//                                        });
//  });
//
//  dirichletValues.fixDOFs([&](auto& basis, auto&& dirichletFlags) {
//    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basis,Dune::Indices::_0, 1),
//                                        [&](auto&& localIndex, auto&& localView, auto&& intersection) {
//                                          if (std::abs(intersection.geometry().center()[0]) > 10 - 1e-8)
//                                            dirichletFlags[localView.index(localIndex)] = true;
//                                        });
//  });

  auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;

  auto req = Ikarus::FErequirements<std::reference_wrapper<MultiTypeVector>>()
      .addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto residualFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getReducedVector(req);
  };

  auto KFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getReducedMatrix(req);
  };

  auto energyFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getScalar(req);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(functions(residualFunction, KFunction), parameter(x, lambda));

  const double gradTol = 1e-16;

  // auto tr = Ikarus::makeTrustRegion(nonLinOp);
  // tr->setup({.verbosity = 1,
  //            .maxiter   = 1000,
  //            .grad_tol  = gradTol,
  //            .corr_tol  = 1e-8,  // everything should converge to the gradient tolerance
  //            .useRand   = false,
  //            .rho_reg   = 1e8,
  //            .Delta0    = 1});
  auto updateFunction = std::function([&](MultiTypeVector &multiTypeVector, const Eigen::VectorXd &d_) {
    auto dFull = sparseAssembler.createFullVector(d_);
//    std::cout<<"dFull"<<std::endl;
//    std::cout<<dFull<<std::endl;
    multiTypeVector += dFull;
  });
  auto linSolver               = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_UmfPackLU);
  auto tr                      = Ikarus::makeNewtonRaphson(nonLinOp, std::move(linSolver),std::move(updateFunction));

  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  tr->subscribeAll(nonLinearSolverObserver);

  auto basis3D       = Ikarus::makeBasis(gridView, composite(power<3>(scalarMidSurfBasis, BlockedInterleaved()),power<3>(scalaarDirectorBasis, BlockedInterleaved()),
                                                           BlockedLexicographic{}));
  auto blockedmidSurfaceBasis2 = Dune::Functions::subspaceBasis(basis3D.untouched(),Dune::Indices::_0);
  auto blockeddirectorBasis2 = Dune::Functions::subspaceBasis(basis3D.untouched(),Dune::Indices::_1);
  Dune::Vtk::Shell3DDataCollector dataCollector1(gridView,thickness,Dune::RefinementIntervals(plotInPlaneRefine));
  auto f=[thickness](auto&& m, auto&& t,auto&& t0, auto&& x)
  {
    //    std::cout<<"m"<<std::endl;
    //    std::cout<<m<<std::endl;
    //    std::cout<<"t"<<std::endl;
    //    std::cout<<t<<std::endl;
    //    std::cout<<"t0"<<std::endl;
    //    std::cout<<t0<<std::endl;
    //    std::cout<<"m+(2*x-1)*(t-t0)* thickness / 2.0"<<std::endl;
    //    std::cout<<m+(2*x-1)*(t-t0)* thickness / 2.0<<std::endl;
    return m+(2*x-1)*(t-t0)* thickness / 2.0;
  };
  auto resReq = Ikarus::ResultRequirements< Ikarus::FErequirements<std::reference_wrapper<MultiTypeVector>>>()
      .insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, x)
      .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
      .addResultRequest(ResultType::cauchyStress);
  auto director0 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(blockeddirectorBasis2,  x0);


  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis)>,MultiTypeVector>>(
      basis, x,[&](auto& writer, auto& basis,auto& xL, auto& prefixString, int step){

        auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(blockedmidSurfaceBasis2,  x);
        auto director = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(blockeddirectorBasis2,  x);

        writer.addVertexData(disp, {"displacements", Dune::VTK::FieldInfo::Type::vector, 3});
        writer.addVertexData(director, {"director", Dune::VTK::FieldInfo::Type::vector, 3});
        writer.write(prefixString+ std::to_string(step));
        Dune::VtkUnstructuredGridWriterMod writer2(dataCollector1, Dune::Vtk::FormatTypes::ASCII);


        auto compf = Dune::Functions::ComposedGridFunctionMod(f,disp,director,director0);
        //        std::cout<<"T3"<<std::endl;
        auto resultFunction = std::make_shared<ResultFunction3D<ElementType>>(&fes, resReq);
//        auto resultFunction2 = std::make_shared<Dune::VTKFunctionMod<GridView>const>(resultFunction);

        writer2.addPointData(compf, Dune::Vtk::FieldInfo{"displacements", Dune::Vtk::RangeTypes::VECTOR, 3});
        writer2.addPointData(Dune::Vtk::FunctionMod<GridView>(resultFunction));
        //        std::cout<<"T4"<<std::endl;
        const std::string name3d= "PureBending3D"+std::to_string(step);
        std::cout<<name3d<<std::endl;
        writer2.write(name3d);
//        std::cout<<"T5"<<std::endl;
      },
      2);
  vtkWriter->setFileNamePrefix("PureBending");
//  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 3);

  auto lc = Ikarus::LoadControl(tr, loadSteps, {0, 1});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run();







  std::cout << std::setprecision(16) << std::ranges::max(d) << std::endl;
  t.check(Dune::FloatCmp::eq(0.2957393081676369, std::ranges::max(d)))
      << std::setprecision(16) << "The maximum displacement is " << std::ranges::max(d);
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
//  checkFEByAutoDiff<RMSHELLSB>("RMSHELLSB");

//  checkFEByAutoDiff<KLSHELLSB>("KLSHELLSB");
  NonLinearElasticityLoadControlNRandTRforRMShell();
}
