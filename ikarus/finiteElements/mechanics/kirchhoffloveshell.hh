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

template<typename Basis_,  typename FERequirements_ = FErequirements<>, bool useEigenRef_ = false>
class KirchhoffLoveShell : public PowerBasisFE<typename Basis_::FlatBasis> {
 public:
  using Basis = Basis_;
  using FlatBasis = typename Basis::FlatBasis;
  using BasePowerFE = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
  using FERequirementType = FERequirements_;
  using ResultRequirementsType = ResultRequirements<FERequirementType>;
  using LocalView = typename FlatBasis::LocalView;
  using Element = typename LocalView::Element;
  using Geometry = typename Element::Geometry;
  using GridView = typename FlatBasis::GridView;
  static constexpr int useEigenRef = useEigenRef_;
  using Traits = TraitsFromLocalView<LocalView, useEigenRef>;
  static constexpr int myDim = Traits::mydim;
  static constexpr int worlddim = Traits::worlddim;

  template<typename VolumeLoad = LoadDefault, typename NeumannBoundaryLoad = LoadDefault>
  KirchhoffLoveShell(const Basis &globalBasis, const typename LocalView::Element &element,
                     FESettings p_feSettings, VolumeLoad p_volumeLoad = {},
                     const BoundaryPatch<GridView> *p_neumannBoundary = nullptr,
                     NeumannBoundaryLoad p_neumannBoundaryLoad = {})
      : BasePowerFE(globalBasis.flat(), element), neumannBoundary{p_neumannBoundary}, fESettings{std::move(p_feSettings)} {
    const auto &simulationFlag = fESettings.request<int>("simulationFlag");
    if (simulationFlag==0)
      membraneStrain = DefaultMembraneStrain();
    else if (simulationFlag==1)
      membraneStrain = CASMembraneStrain();

    this->localView().bind(element);
    auto &first_child = this->localView().tree().child(0);
    const auto &fe = first_child.finiteElement();
    numberOfNodes = fe.size();
    dispAtNodes.resize(fe.size());
    dispAtNodesDual.resize(fe.size());
    dispAtNodesDual2nd.resize(fe.size());
    order = 2*(fe.localBasis().order());
    localBasis = Dune::CachedLocalBasis(fe.localBasis());
    if constexpr (requires { this->localView().element().impl().getQuadratureRule(order); })
      if (this->localView().element().impl().isTrimmed())
        localBasis.bind(this->localView().element().impl().getQuadratureRule(order), Dune::bindDerivatives(0, 1, 2));
      else
        localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                        Dune::bindDerivatives(0, 1, 2));
    else
      localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                      Dune::bindDerivatives(0, 1, 2));

    if constexpr (!std::is_same_v<VolumeLoad, LoadDefault>) volumeLoad = p_volumeLoad;
    if constexpr (!std::is_same_v<NeumannBoundaryLoad, LoadDefault>) neumannBoundaryLoad = p_neumannBoundaryLoad;

    assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
               && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
  }

 public:
  template<typename ScalarType>
  auto getDisplacementFunction(const FERequirementType &par,
                               const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
    const auto &d = par.getGlobalSolution(FESolutions::displacement);
    auto geo = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
        this->localView().element().geometry());
    Dune::BlockVector<Dune::RealTuple<ScalarType, Traits::worlddim>> disp(dispAtNodes.size());

    if (dx)
      for (auto i = 0U; i < disp.size(); ++i)
        for (auto k2 = 0U; k2 < worlddim; ++k2)
          disp[i][k2]
              = dx.value()[i*worlddim + k2]
              + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
    else
      for (auto i = 0U; i < disp.size(); ++i)
        for (auto k2 = 0U; k2 < worlddim; ++k2)
          disp[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];

    Dune::StandardLocalFunction uFunction(localBasis, disp, geo);

    return uFunction;
  }

  inline double calculateScalar(const FERequirementType &par) const { return calculateScalarImpl<double>(par); }

  void calculateAt(const ResultRequirementsType &req, const Dune::FieldVector<double, Traits::mydim> &local,
                   ResultTypeMap<double> &result) const {
    DUNE_THROW(Dune::NotImplemented, "No results are implemented");
  }

  inline void calculateVector(const FERequirementType &par, typename Traits::template VectorType<> force) const {
    calculateVectorImpl<double>(par, force);
  }

  inline void calculateMatrix(const FERequirementType &par, typename Traits::template MatrixType<> K) const {
    calculateMatrixImpl<double>(par, K);
  }

  Dune::CachedLocalBasis<
      std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
      localBasis;
  std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim> &,
                                                        const double &)>
      volumeLoad;
  std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim> &,
                                                        const double &)>
      neumannBoundaryLoad;
  const BoundaryPatch<GridView> *neumannBoundary;
  mutable Dune::BlockVector<Dune::RealTuple<double, Traits::worlddim>> dispAtNodes;
  mutable Dune::BlockVector<Dune::RealTuple<autodiff::dual, Traits::worlddim>> dispAtNodesDual;
  mutable Dune::BlockVector<Dune::RealTuple<autodiff::dual2nd, Traits::worlddim>> dispAtNodesDual2nd;
  FESettings fESettings;
  size_t numberOfNodes{0};
  int order{};
  mutable MembraneStrain membraneStrain;

  auto calc3DMetric(const auto& J,const auto& H,double zeta) const
  {
    const auto A1 = J.row(0).transpose();
    const auto A2 = J.row(1).transpose();
    const auto A1xA2 = (A1.cross(A2)).eval();
    const auto A3norm = A1xA2.norm();
    const auto A3 = (A1xA2/A3norm).eval();
    const auto dA1_dT1 = H.row(0).transpose();
    const auto dA2_dT2 = H.row(1).transpose();
    const auto dA2_dT1 = H.row(2).transpose();
    const auto dA1_dT2 = H.row(2).transpose();

    const auto  d_A1crossA2_dT1 = (dA1_dT1.cross(A2)+ A1.cross(dA2_dT1)).eval();
    const auto  d_A1crossA2_dT2 = (dA1_dT2.cross(A2)+ A1.cross(dA2_dT2)).eval();

    const Eigen::Matrix<double, 3, 3>
        P = 1.0/A3norm*(Eigen::Matrix<double, 3, 3>::Identity() - A3*A3.transpose());

    const auto  dA3_dT1 = (P*d_A1crossA2_dT1).eval();
    const auto  dA3_dT2 = (P*d_A1crossA2_dT2).eval();

    const auto G1 = (A1 + zeta*dA3_dT1).eval();
    const auto G2 = (A2 + zeta*dA3_dT2).eval();
    const auto G3 = A3;

    Eigen::Matrix<double,3,3> J3D;
    J3D.col(0)=G1;
    J3D.col(1)=G2;
    J3D.col(2)=G3;

    Eigen::Matrix<double,3,3> G= J3D.transpose()*J3D; //here quadratic parts of zeta appear!
    return G;
  }

  friend StressBasedShell<KirchhoffLoveShell>;
 protected:
  mutable std::vector<Dune::FieldVector<double, 2>> lagrangePoints;
  mutable std::vector<double> out;

  auto computeMaterialAndStrains(const Dune::FieldVector<double, 2> &gpPos,
                                 int gpIndex,
                                 const Geometry &geo,
                                 const auto &uFunction) const {
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    const auto [X, Jd, Hd] = geo.impl().zeroFirstAndSecondDerivativeOfPosition(gpPos);
    const auto J = toEigen(Jd);
    const auto H = toEigen(Hd);
    const Eigen::Matrix<double, 2, 2> A = J*J.transpose();
    Eigen::Matrix<double, 3, 3> G;
    G.setZero();
    G.block<2, 2>(0, 0) = A;
    G(2, 2) = 1;
    const Eigen::Matrix<double, 3, 3> GInv = G.inverse();
    const auto C = materialTangent(GInv);

    Eigen::Vector3<ScalarType> epsV=membraneStrain.value(gpPos,geo,uFunction);

    const auto &Ndd = localBasis.evaluateSecondDerivatives(gpIndex);
    const auto uasMatrix = Dune::viewAsEigenMatrixAsDynFixed(uFunction.coefficientsRef());
    const auto hessianu = Ndd.transpose().template cast<ScalarType>()*uasMatrix;
    const Eigen::Matrix3<ScalarType> h = H + hessianu;
    const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
        uFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll), Dune::on(Dune::DerivativeDirections::referenceElement)));
    const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();
    const Eigen::Vector3<ScalarType> a3N = (j.row(0).cross(j.row(1)));
    const Eigen::Vector3<ScalarType> a3 = a3N.normalized();
    Eigen::Vector<ScalarType, 3> bV = h*a3;
    bV(2) *= 2;  // Voigt notation requires the two here
    const auto BV = toVoigt(toEigen(geo.impl().secondFundamentalForm(gpPos)));
    const auto kappaV = (BV - bV).eval();
    return std::make_tuple(C, epsV, kappaV, j, J, h,H, a3N, a3);
  }

  template<typename ScalarType>
  auto calculateScalarImpl(const FERequirementType &par, const std::optional<const Eigen::VectorX<ScalarType>> &dx
  = std::nullopt) const -> ScalarType {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = getDisplacementFunction(par, dx);
    const auto &lambda = par.getParameter(FEParameter::loadfactor);
    const auto geo = this->localView().element().geometry();
    ScalarType energy = 0.0;
    this->membraneStrain.pre(geo,uFunction);

    const auto &thickness_ = fESettings.request<double>("thickness");
    for (const auto &[gpIndex, gp]: uFunction.viewOverIntegrationPoints()) {

      const auto
          [C, epsV, kappaV, j, J, h,H, a3N, a3] = computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);

      const ScalarType membraneEnergy = 0.5*thickness_*epsV.dot(C*epsV);
      const ScalarType bendingEnergy = 0.5*Dune::power(thickness_, 3)/12.0*kappaV.dot(C*kappaV);
      energy += (membraneEnergy + bendingEnergy)*geo.integrationElement(gp.position())*gp.weight();
    }

      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
          const auto u                                       = uFunction.evaluate(gpIndex);
          const Eigen::Vector<double, Traits::worlddim> fExt = volumeLoad(toEigen(geo.global(gp.position())), lambda);
          energy -= u.dot(fExt) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e., neumann boundary
      if (not neumannBoundary and not neumannBoundaryLoad) return energy;
      forEachInterSectionIntegrationPoint(this->localView().element(),neumannBoundary,order,
                                                  [&](auto& quadPos,auto&& globalPos,auto& intElement){
                                                                          const auto neumannValue
                                                                              = neumannBoundaryLoad(globalPos, lambda);
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
    const auto uFunction = getDisplacementFunction(par, dx);
    const auto &lambda = par.getParameter(FEParameter::loadfactor);
    const auto geo = this->localView().element().geometry();
    this->membraneStrain.pre(geo,uFunction);

    const auto &thickness_ = fESettings.request<double>("thickness");

    // Internal forces
    for (const auto &[gpIndex, gp]: uFunction.viewOverIntegrationPoints()) {

      const auto
          [C, epsV, kappaV, jE, J, h,H, a3N, a3] = computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);
      const Eigen::Vector<ScalarType, 3> membraneForces = thickness_*C*epsV;
      const Eigen::Vector<ScalarType, 3> moments = Dune::power(thickness_, 3)/12.0*C*kappaV;

      const auto &Nd = localBasis.evaluateJacobian(gpIndex);
      const auto &Ndd = localBasis.evaluateSecondDerivatives(gpIndex);
      for (size_t i = 0; i < numberOfNodes; ++i) {
        Eigen::Matrix<ScalarType, 3, 3> bopIMembrane = membraneStrain.derivative(gp.position(),jE, Nd, geo,uFunction,localBasis, i);
        Eigen::Matrix<ScalarType, 3, 3> bopIBending = bopBending(jE, h, Nd, Ndd, i, a3N, a3);
        force.template segment<3>(3*i) +=
            bopIMembrane.transpose()*membraneForces*geo.integrationElement(gp.position())*gp.weight();
        force.template segment<3>(3*i) +=
            bopIBending.transpose()*moments*geo.integrationElement(gp.position())*gp.weight();
      }
    }

     //External forces volume forces over the domain
          if (volumeLoad) {
            const auto u = getDisplacementFunction(par, dx);
            for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
              Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigen(geo.global(gp.position())), lambda);
              for (size_t i = 0; i < numberOfNodes; ++i) {
                const auto udCi = uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
                force.template segment<worlddim>(worlddim * i)
                    -= udCi * fext * geo.integrationElement(gp.position()) * gp.weight();
              }
            }
          }

          // External forces, boundary forces, i.e. at the Neumann boundary
          if (not neumannBoundary) return;
          forEachInterSectionIntegrationPoint(this->localView().element(),neumannBoundary,order,
                                      [&](auto& quadPos,auto&& globalPos,auto& intElement){
                                        for (size_t i = 0; i < numberOfNodes; ++i) {
                                          const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

                                          // Value of the Neumann data at the current position
                                          auto neumannValue
                                              = neumannBoundaryLoad(globalPos, lambda);
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
    const auto uFunction = getDisplacementFunction(par, dx);
    const auto &lambda = par.getParameter(FEParameter::loadfactor);
    const auto geo = this->localView().element().geometry();
    this->membraneStrain.pre(geo,uFunction);

    const auto &thickness_ = fESettings.request<double>("thickness");

    // Internal forces
    for (const auto &[gpIndex, gp]: uFunction.viewOverIntegrationPoints()) {
      const auto intElement = geo.integrationElement(gp.position())*gp.weight();
      const auto
          [C, epsV, kappaV, jE, J, h,H, a3N, a3] = computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);
      const Eigen::Vector<ScalarType, 3> membraneForces = thickness_*C*epsV;
      const Eigen::Vector<ScalarType, 3> moments = Dune::power(thickness_, 3)/12.0*C*kappaV;

      const auto &Nd = localBasis.evaluateJacobian(gpIndex);
      const auto &Ndd = localBasis.evaluateSecondDerivatives(gpIndex);
      for (size_t i = 0; i < numberOfNodes; ++i) {
        Eigen::Matrix<ScalarType, 3, 3> bopIMembrane = membraneStrain.derivative(gp.position(),jE, Nd,geo,uFunction,localBasis, i);
        Eigen::Matrix<ScalarType, 3, 3> bopIBending = bopBending(jE, h, Nd, Ndd, i, a3N, a3);
        for (size_t j = i; j < numberOfNodes; ++j) {
          Eigen::Matrix<ScalarType, 3, 3> bopJMembrane = membraneStrain.derivative(gp.position(),jE, Nd,geo,uFunction,localBasis, j);
          Eigen::Matrix<ScalarType, 3, 3> bopJBending = bopBending(jE, h, Nd, Ndd, j, a3N, a3);
          K.template block<3, 3>(3*i, 3*j) += thickness_*bopIMembrane.transpose()*C*bopJMembrane*intElement;
          K.template block<3, 3>(3*i, 3*j) +=
              Dune::power(thickness_, 3)/12.0*bopIBending.transpose()*C*bopJBending*intElement;

          Eigen::Matrix<ScalarType, 3, 3> kgMembraneIJ = membraneStrain.secondDerivative(gp.position(),Nd,geo,uFunction,localBasis, membraneForces, i, j);
          Eigen::Matrix<ScalarType, 3, 3> kgBendingIJ = kgBending(jE, h, Nd, Ndd, a3N, a3, moments, i, j);
          K.template block<3, 3>(3*i, 3*j) += kgMembraneIJ*intElement;
          K.template block<3, 3>(3*i, 3*j) += kgBendingIJ*intElement;
        }
      }
    }
    K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
  }

 private:
  template<typename ScalarType>
  Eigen::Matrix<ScalarType, 3, 3> bopMembrane(const Eigen::Matrix<ScalarType, 2, 3> &Jcur, const auto &dN,
                                              const int node) const {
    Eigen::Matrix<ScalarType, 3, 3> bop;
    bop.row(0) = Jcur.row(0)*dN(node, 0);
    bop.row(1) = Jcur.row(1)*dN(node, 1);
    bop.row(2) = Jcur.row(0)*dN(node, 1) + Jcur.row(1)*dN(node, 0);

    return bop;
  }

  template<typename ScalarType>
  Eigen::Matrix<ScalarType, 3, 3> kgMembrane(const auto &dN,
                                             const Eigen::Vector3<ScalarType> &S, int I, int J) const {

    const auto &dN1i = dN(I, 0);
    const auto &dN1j = dN(J, 0);
    const auto &dN2i = dN(I, 1);
    const auto &dN2j = dN(J, 1);
    const ScalarType NS = dN1i*dN1j*S[0] + dN2i*dN2j*S[1] + (dN1i*dN2j + dN2i*dN1j)*S[2];
    Eigen::Matrix<ScalarType, 3, 3> kg = Eigen::Matrix<double, 3, 3>::Identity()*NS;
    return kg;
  }

  template<typename ScalarType>
  Eigen::Matrix<ScalarType, 3, 3> kgBending(const Eigen::Matrix<ScalarType, 2, 3> &jcur,
                                            const Eigen::Matrix3<ScalarType> &h,
                                            const auto &dN,
                                            const auto &ddN,
                                            const Eigen::Vector3<ScalarType> &a3N,
                                            const Eigen::Vector3<ScalarType> &a3,
                                            const Eigen::Vector3<ScalarType> &S,
                                            int I,
                                            int J) const {
    Eigen::Matrix<ScalarType, 3, 3> kg;
    kg.setZero();

    const auto &dN1i = dN(I, 0);
    const auto &dN1j = dN(J, 0);
    const auto &dN2i = dN(I, 1);
    const auto &dN2j = dN(J, 1);

    const Eigen::Matrix<ScalarType, 3, 3>
        P = 1.0/a3N.norm()*(Eigen::Matrix<double, 3, 3>::Identity() - a3*a3.transpose());

    const auto a1dxI =
        Eigen::Matrix<double, 3, 3>::Identity()*dN1i; // the auto here allows the exploitation of the identity matrices,
    // due to Eigen's expression templates
    const auto a2dxI = Eigen::Matrix<double, 3, 3>::Identity()*dN2i;
    const auto a1dxJ = Eigen::Matrix<double, 3, 3>::Identity()*dN1j;
    const auto a2dxJ = Eigen::Matrix<double, 3, 3>::Identity()*dN2j;
    const auto a1 = jcur.row(0);
    const auto a2 = jcur.row(1);
    const Eigen::Matrix<ScalarType, 3, 3> a3NdI = a1dxI.colwise().cross(a2) - a2dxI.colwise().cross(a1);
    const Eigen::Matrix<ScalarType, 3, 3> a3NdJ = a1dxJ.colwise().cross(a2) - a2dxJ.colwise().cross(a1);
    Eigen::Matrix<ScalarType, 3, 3> a3dI = P*a3NdI;
    Eigen::Matrix<ScalarType, 3, 3> a3dJ = P*a3NdJ;
    for (int i = 0; i < 3; ++i) {

      const auto a_albe = h.row(i).transpose();
      const auto &ddNI = ddN(I, i);
      const auto &ddNJ = ddN(J, i);
      Eigen::Vector3<ScalarType> vecd = P*a_albe;

      Eigen::Matrix<ScalarType, 3, 3> a3Ndd = 1.0/a3N.squaredNorm()
          *((3*a3*a3.transpose() - Eigen::Matrix<double, 3, 3>::Identity())*(a3.dot(a_albe))
              - a_albe*a3.transpose() - a3*a_albe.transpose());

      Eigen::Matrix<ScalarType, 3, 3>
          secondDerivativeDirectorIJ = skew(((dN2i*dN1j - dN1i*dN2j)*vecd).eval());
      kg -= (a3NdI.transpose()*a3Ndd*a3NdJ + secondDerivativeDirectorIJ + (ddNI*a3dJ + ddNJ*a3dI.transpose()))*S[i]
          *(i==2 ? 2 : 1);
    }

//      const ScalarType NS = dN1i * dN1j * S[0] + dN2i * dN2j * S[1] + (dN1i * dN2j + dN2i * dN1j) * S[2];
//      Eigen::Matrix<ScalarType, 3, 3> kg= Eigen::Matrix<double, 3, 3>::Identity()*NS;
    return kg;
  }

  template<typename ScalarType>
  Eigen::Matrix<ScalarType, 3, 3> bopBending(const Eigen::Matrix<ScalarType, 2, 3> &jcur,
                                             const Eigen::Matrix3<ScalarType> &h,
                                             const auto &dN,
                                             const auto &ddN,
                                             const int node,
                                             const Eigen::Vector3<ScalarType> &a3N,
                                             const Eigen::Vector3<ScalarType> &a3) const {
    const Eigen::Matrix<ScalarType, 3, 3>
        a1dxI = Eigen::Matrix<double, 3, 3>::Identity()*dN(node, 0); // this should be double
    // but the cross-product below complains otherwise
    const Eigen::Matrix<ScalarType, 3, 3> a2dxI = Eigen::Matrix<double, 3, 3>::Identity()*dN(node, 1);
    const auto a1 = jcur.row(0);
    const auto a2 = jcur.row(1);
    const Eigen::Matrix<ScalarType, 3, 3>
        a3NdI = a1dxI.colwise().cross(a2) - a2dxI.colwise().cross(a1); // minus needed since order has
    // to be swapped to get column-wise cross product working
    const Eigen::Matrix<ScalarType, 3, 3>
        a3d1 = 1.0/a3N.norm()*(Eigen::Matrix<double, 3, 3>::Identity() - a3*a3.transpose())*a3NdI;

    Eigen::Matrix<ScalarType, 3, 3> bop = -(h*a3d1 + (a3*ddN.row(node)).transpose());
    bop.row(2) *= 2;

    return bop;
  }

  Eigen::Matrix<double, 3, 3> materialTangent(const Eigen::Matrix<double, 3, 3> &Aconv) const {
    const auto &emod_ = fESettings.request<double>("youngs_modulus");
    const auto &nu_ = fESettings.request<double>("poissons_ratio");
    const double lambda = emod_*nu_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    const double mu = emod_/(2.0*(1.0 + nu_));
    const double lambdbar = 2.0*lambda*mu/(lambda + 2.0*mu);
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>> moduli;
    const auto AconvT = TensorCast(Aconv, std::array<Eigen::Index, 2>({3, 3}));
    moduli = lambdbar*dyadic(AconvT, AconvT).eval() + 2.0*mu*symmetricFourthOrder<double>(Aconv, Aconv);

    auto C = toVoigt(moduli);
    Eigen::Matrix<double, 3, 3> C33 = C({0, 1, 5}, {0, 1, 5});
    return C33;
  }
};


}  // namespace Ikarus
