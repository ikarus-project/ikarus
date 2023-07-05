// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "externalload.hh"
#include "fesettings.hh"

#include <utility>

#include <dune/common/classname.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/localfefunctions/manifolds/realTuple.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feBases/powerBasisFE.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/finiteElements/mechanics/membranestrains.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/utils/tensorUtils.hh>
#include <ikarus/utils/tensorproductquadrule.hh>
namespace Ikarus {


  template<typename ResultantBasedShell>
  struct StressBasedShell : private ResultantBasedShell
  {
    Dune::QuadratureRule<double,3> rule;
    template <typename... Args>
    explicit StressBasedShell(Args&&... args) : ResultantBasedShell{std::forward<Args>(args)...} {

      const auto& twoDRule = Dune::QuadratureRules<double,2>::rule(Dune::GeometryTypes::quadrilateral,this->order);
      const auto& oneDRule = Dune::QuadratureRules<double,1>::rule(Dune::GeometryTypes::line,2);

      numberOfThicknessIntegrationPoints = oneDRule.size();
      rule = Ikarus::tensorProductQuadrature(twoDRule,oneDRule);
//      fESettings = ResultantBasedShell::getfFESettings();
//      order = ResultantBasedShell::getOrder();
//      numberOfNodes = ResultantBasedShell::getNumberOfNodes();
//      membraneStrain = ResultantBasedShell::getMembraneStrain();
//      std::cout<<"Constructor: "<<rule.size()<<std::endl;

    }

    int numberOfThicknessIntegrationPoints;
//    FESettings fESettings;
//    size_t numberOfNodes{0};
//    int order{};
//    mutable MembraneStrain membraneStrain;

    using GlobalIndex = typename ResultantBasedShell::GlobalIndex;
    using ResultantBasedShell::size ;
    using ResultantBasedShell::globalFlatIndices ;
    using ResultantBasedShell::localView ;
    using Basis = typename ResultantBasedShell::Basis;
    using FlatBasis = typename Basis::FlatBasis;
    using BasePowerFE = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType = typename ResultantBasedShell::FERequirementType;
    using ResultRequirementsType = ResultRequirements<FERequirementType>;
    using LocalView = typename FlatBasis::LocalView;
    using Element = typename LocalView::Element;
    using Geometry = typename Element::Geometry;
    using GridView = typename FlatBasis::GridView;
    static constexpr int useEigenRef = ResultantBasedShell::useEigenRef;
    using Traits = TraitsFromLocalView<LocalView, useEigenRef>;
    static constexpr int myDim = Traits::mydim;
    static constexpr int worlddim = Traits::worlddim;

    inline void calculateMatrix(const FERequirementType &par, typename Traits::template MatrixType<> K) const {
      calculateMatrixImpl<double>(par, K);
    }

    inline void calculateVector(const FERequirementType &par, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(par, force);
    }
    inline double calculateScalar(const FERequirementType &par) const { return calculateScalarImpl<double>(par); }
    template<typename ScalarType>
    auto calculateScalarImpl(const FERequirementType &par, const std::optional<const Eigen::VectorX<ScalarType>> &dx
                                                           = std::nullopt) const -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = this->getDisplacementFunction(par, dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      const auto geo = this->localView().element().geometry();
      ScalarType energy = 0.0;
      this->membraneStrain.pre(geo,uFunction);

      const auto &thickness_ = this->fESettings.template request<double>("thickness");
      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;
        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};
        const auto
            [C, epsV, kappaV, j, J, h,H, a3N, a3] = this->computeMaterialAndStrains(gp2DPos, gpIndex2D, geo, uFunction);

        const auto G = this->calc3DMetric(J,H,zeta);
//        std::cout<<"G: "<<G<<std::endl;
        const auto Ginv = G.inverse().eval();

        const auto C3D = this->materialTangent(Ginv);
//        std::cout<<"C3D: "<<C3D<<std::endl;
        const auto strainsV= (epsV+ zeta*kappaV).eval();
//        std::cout<<"strainsV: "<<strainsV<<std::endl;
        const ScalarType energyVal = 0.5*strainsV.dot(C3D*strainsV);
        energy += (energyVal)*geo.integrationElement(gp2DPos)*gp.weight()*2*thickness_/2.0;
        ++gpIndex;
      }

      if (this->volumeLoad) {
        for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
          const auto u                                       = uFunction.evaluate(gpIndex);
          const Eigen::Vector<double, Traits::worlddim> fExt = this->volumeLoad(toEigen(geo.global(gp.position())), lambda);
          energy -= u.dot(fExt) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e., neumann boundary
      if (not this->neumannBoundary and not this->neumannBoundaryLoad) return energy;
      forEachInterSectionIntegrationPoint(this->localView().element(),this->neumannBoundary,this->order,
                                          [&](auto& quadPos,auto&& globalPos,auto& intElement){
                                            const auto neumannValue
                                                = this->neumannBoundaryLoad(globalPos, lambda);
                                            const auto u = uFunction.evaluate(quadPos);
                                            energy -= neumannValue.dot(u) * intElement;
                                          });
      return energy;
    }

    template<typename ScalarType>
    void calculateVectorImpl(const FERequirementType &par, typename Traits::template VectorType<ScalarType> force,
                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
      force.setZero();
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = this->getDisplacementFunction(par, dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      const auto geo = this->localView().element().geometry();

      const auto &thickness_ = this->fESettings.template request<double>("thickness");
      this->membraneStrain.pre(geo,uFunction);

      // Internal forces
      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;
        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};
        const auto
            [C, epsV, kappaV, j, J, h,H, a3N, a3] = this->computeMaterialAndStrains(gp2DPos, gpIndex2D, geo, uFunction);
        const auto G = this->calc3DMetric(J,H,zeta);
        //        std::cout<<"G: "<<G<<std::endl;
        const auto Ginv = G.inverse().eval();

        const auto C3D = this->materialTangent(Ginv);
        const auto strainsV= (epsV+ zeta*kappaV).eval();
        const auto S = (C3D*strainsV).eval();

        const auto &Nd = this->localBasis.evaluateJacobian(gpIndex2D);
        const auto &Ndd = this->localBasis.evaluateSecondDerivatives(gpIndex2D);
        for (size_t i = 0; i < this->numberOfNodes; ++i) {
          Eigen::Matrix<ScalarType, 3, 3> bopIMembrane = this->membraneStrain.derivative(gp2DPos,j, Nd,geo,uFunction,this->localBasis, i);
          Eigen::Matrix<ScalarType, 3, 3> bopIBending = this->bopBending(j, h, Nd, Ndd, i, a3N, a3);
          Eigen::Matrix<ScalarType, 3, 3> bopI = bopIMembrane+zeta*bopIBending;
          force.template segment<3>(3*i) +=
              bopI.transpose()*S*geo.integrationElement(gp2DPos)*gp.weight()*2*thickness_/2.0;
          // the first two fixes the change of the integration mapping from 0..1 to -1..1,
          // and the h/2 factor is the factor for the correct thickness
        }
        ++gpIndex;
      }

      //External forces volume forces over the domain
      if (this->volumeLoad) {
        const auto u = this->getDisplacementFunction(par, dx);
        for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
          Eigen::Vector<double, Traits::worlddim> fext = this->volumeLoad(toEigen(geo.global(gp.position())), lambda);
          for (size_t i = 0; i < this->numberOfNodes; ++i) {
            const auto udCi = uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<worlddim>(worlddim * i)
                -= udCi * fext * geo.integrationElement(gp.position()) * gp.weight();
          }
        }
      }

      // External forces, boundary forces, i.e., at the Neumann boundary
      if (not this->neumannBoundary) return;
      forEachInterSectionIntegrationPoint(this->localView().element(),this->neumannBoundary,this->order,
                                          [&](auto& quadPos,auto&& globalPos,auto& intElement){
                                            for (size_t i = 0; i < this->numberOfNodes; ++i) {
                                              const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

                                              // Value of the Neumann data at the current position
                                              auto neumannValue
                                                  = this->neumannBoundaryLoad(globalPos, lambda);
                                              force.template segment<worlddim>(worlddim * i) -= udCi * neumannValue * intElement;
                                            }
                                          });
    }

    template<typename ScalarType>
    void calculateMatrixImpl(const FERequirementType &par, typename Traits::template MatrixType<ScalarType> K,
                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
      K.setZero();
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = this->getDisplacementFunction(par, dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      const auto geo = this->localView().element().geometry();
      this->membraneStrain.pre(geo,uFunction);

      const auto &thickness_ = this->fESettings.template request<double>("thickness");

      // Internal forces
      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;
        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};
        const auto
          [C, epsV, kappaV, jE, J, h,H, a3N, a3] = this->computeMaterialAndStrains(gp2DPos, gpIndex2D, geo, uFunction);
        const auto G = this->calc3DMetric(J,H,zeta);
        //        std::cout<<"G: "<<G<<std::endl;
        const auto Ginv = G.inverse().eval();

        const auto C3D = this->materialTangent(Ginv);
        const auto strainsV= (epsV+ zeta*kappaV).eval();
        const auto S = (C3D*strainsV).eval();
        const auto intElement =geo.integrationElement(gp2DPos) * gp.weight() * 2 * thickness_ / 2.0;

        const auto &Nd = this->localBasis.evaluateJacobian(gpIndex2D);
        const auto &Ndd = this->localBasis.evaluateSecondDerivatives(gpIndex2D);
        for (size_t i = 0; i < this->numberOfNodes; ++i) {
          Eigen::Matrix<ScalarType, 3, 3> bopIMembrane = this->membraneStrain.derivative(gp2DPos,jE, Nd,geo,uFunction,this->localBasis, i);
          Eigen::Matrix<ScalarType, 3, 3> bopIBending = this->bopBending(jE, h, Nd, Ndd, i, a3N, a3);
          Eigen::Matrix<ScalarType, 3, 3> bopI = bopIMembrane + zeta * bopIBending;
          for (size_t j = i; j < this->numberOfNodes; ++j) {
            Eigen::Matrix<ScalarType, 3, 3> bopJMembrane = this->membraneStrain.derivative(gp2DPos,jE, Nd,geo,uFunction,this->localBasis, j);
            Eigen::Matrix<ScalarType, 3, 3> bopJBending = this->bopBending(jE, h, Nd, Ndd, j, a3N, a3);
            Eigen::Matrix<ScalarType, 3, 3> bopJ = bopJMembrane + zeta * bopJBending;
            Eigen::Matrix<ScalarType, 3, 3> kgMembraneIJ = this->membraneStrain.secondDerivative(gp2DPos,Nd,geo,uFunction,this->localBasis, S, i, j);
            Eigen::Matrix<ScalarType, 3, 3> kgBendingIJ = this->kgBending(jE, h, Nd, Ndd, a3N, a3, S, i, j);
            K.template block<3, 3>(3*i, 3*j) += (bopI.transpose()*C*bopJ+kgMembraneIJ+zeta*kgBendingIJ)*intElement;

          }
        }
        ++gpIndex;
      }
      K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
          }
  };

}  // namespace Ikarus