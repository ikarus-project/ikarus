// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/classname.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/expressions/greenLagrangeStrains.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feBases/powerBasisFE.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/finiteElements/mechanics/materials.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  auto toLocalCartesian(const auto& A,const auto& Jcontravariant,const auto& Jloc)
  {

    auto factor= (Jloc*Jcontravariant.transpose()).eval();
    auto Aloc=(factor*A*factor.transpose()).eval();
    return Aloc;
  }

  auto calculateContravariantBaseVectors(const auto&J,const auto& A3)
  {
    const auto det = J.row(0).dot(J.row(1).cross(A3));
    auto Jcont= J;
    Jcont.row(0)= J.row(1).cross(A3)/det;
    Jcont.row(1)= A3.cross(J.row(0))/det;
    return Jcont;
  }

  template <typename Basis_, typename FERequirements_ = FErequirements<>, bool useEigenRef = false>
  class KirchhoffLoveShell : public PowerBasisFE<typename Basis_::FlatBasis> {
  public:
    using Basis                      = Basis_;
    using FlatBasis                  = typename Basis::FlatBasis;
    using BasePowerFE                = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType          = FERequirements_;
    using ResultRequirementsType     = ResultRequirements<FERequirementType>;
    using LocalView                  = typename FlatBasis::LocalView;
    using Element                    = typename LocalView::Element;
    using Geometry                   = typename Element::Geometry;
    using GridView                   = typename FlatBasis::GridView;
    using Traits                     = TraitsFromLocalView<LocalView, useEigenRef>;
    static constexpr int myDim       = Traits::mydim;
    static constexpr int worlddim       = Traits::worlddim;

    template <typename VolumeLoad = LoadDefault, typename NeumannBoundaryLoad = LoadDefault>
    KirchhoffLoveShell(const Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu, double thickness,
                     VolumeLoad p_volumeLoad = {}, const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
                     NeumannBoundaryLoad p_neumannBoundaryLoad = {})
        : BasePowerFE(globalBasis.flat(), element), neumannBoundary{p_neumannBoundary}, emod_{emod}, nu_{nu},thickness_{thickness} {
      this->localView().bind(element);
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes     = fe.size();
      dispAtNodes.resize(fe.size());
      order = 2 * (fe.localBasis().order());
      localBasis      = Dune::CachedLocalBasis(fe.localBasis());
      if constexpr (requires { this->localView().element().impl().getQuadratureRule(order); })
        if (this->localView().element().impl().isTrimmed())
          localBasis.bind(this->localView().element().impl().getQuadratureRule(order), Dune::bindDerivatives(0, 1,2));
        else
          localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                          Dune::bindDerivatives(0, 1,2));
      else
        localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                        Dune::bindDerivatives(0, 1,2));

      if constexpr (!std::is_same_v<VolumeLoad, LoadDefault>) volumeLoad = p_volumeLoad;
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, LoadDefault>) neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");

      Cmat<< 1,nu,0,
          nu,1,0,
          0,0,(1-nu)/2.0;
      Cmat*= emod/(1-nu*nu);
    }

  public:
    template <typename ScalarType = double>
    auto getDisplacementFunction(const FERequirementType& par,
                                 const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);

      Dune::BlockVector<Dune::RealTuple<ScalarType, Traits::worlddim>> disp(dispAtNodes.size());
      if (dx)
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            disp[i][k2] = dx.value()[i * worlddim + k2]
                          + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      else
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            disp[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];

      auto geo = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          this->localView().element().geometry());
      Dune::StandardLocalFunction uFunction(localBasis, disp, geo);
      return std::make_pair(uFunction,disp);
    }


    inline double calculateScalar(const FERequirementType& par) const { return calculateScalarImpl<double>(par); }

    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
//      using namespace Dune::DerivativeDirections;
//      using namespace Dune;
//
//      const auto uFunction = getDisplacementFunction(req.getFERequirements());
//      const auto H         = uFunction.evaluateDerivative(local, Dune::wrt(spatialAll), Dune::on(gridElement));
//      const auto E         = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
//      const auto EVoigt    = toVoigt(E);
//      auto PK2             = mat.template stresses<StrainTags::greenLagrangian>(EVoigt);
//
//      typename ResultTypeMap<double>::ResultArray resultVector;
//      if (req.isResultRequested(ResultType::PK2Stress)) {
//        resultVector.resizeLike(PK2);
//        resultVector = PK2;
//        result.insertOrAssignResult(ResultType::PK2Stress, resultVector);
//      }
    }

    Dune::CachedLocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        volumeLoad;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        neumannBoundaryLoad;
    const BoundaryPatch<GridView>* neumannBoundary;
    mutable Dune::BlockVector<Dune::RealTuple<double, Traits::dimension>> dispAtNodes;
    Eigen::Matrix<double,3,3> Cmat;
    double emod_;
    double nu_;
    double thickness_;
    size_t numberOfNodes{0};
    int order{};

  protected:
    template <typename ScalarType>
    auto calculateScalarImpl(const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                                           = std::nullopt) const -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [uFunction,uNodes] = getDisplacementFunction(par, dx);
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto geo       = this->localView().element().geometry();
      ScalarType energy    = 0.0;

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto [X,Jd,Hd] = geo.impl().zeroFirstAndSecondDerivativeOfPosition(gp.position());
        const auto J = toEigen(Jd);
        const auto H = toEigen(Hd);
        const auto A1d1 = H.row(0);
        const auto A2d2 =  H.row(1);
        const auto A1d2 = H.row(2);
        const auto A1 = J.row(0);
        const auto A2 = J.row(1);
        const Eigen::Vector3<double> A3 = (A1.cross(A2)).normalized();
        const auto B2 = toEigen(geo.impl().secondFundamentalForm(gp.position()));
        const Eigen::Matrix<double,2,2> A = J*J.transpose();
        const Eigen::Matrix<ScalarType,3,2> gradu = toEigen(uFunction.evaluateDerivative(gpIndex, wrt(spatialAll,Dune::on(DerivativeDirections::referenceElement))));
        const auto ud1 = gradu.col(0);
        const auto ud2 = gradu.col(1);
        const Eigen::Vector3<ScalarType> a1 = A1.transpose()+ud1;
        const Eigen::Vector3<ScalarType> a2 = A2.transpose()+ud2;
        const Eigen::Vector3<ScalarType> a3 = (a1.cross(a2)).normalized();

        const auto& Ndd= localBasis.evaluateSecondDerivatives(gpIndex);
        const auto cps = geo.impl().controlPoints().directGetAll();

//        const h = H +
        Eigen::Vector3<ScalarType> ad1d1=A1d1;
        Eigen::Vector3<ScalarType> ad2d2=A2d2;
        Eigen::Vector3<ScalarType> ad1d2=A1d2;
        Eigen::Vector3<double> Ad1d1;
        Eigen::Vector3<double> Ad2d2;
        Eigen::Vector3<double> Ad1d2;
        Ad1d1.setZero();
        Ad2d2.setZero();
        Ad1d2.setZero();
        for (int i = 0; i < numberOfNodes; ++i) {
          ad1d1+=Ndd.col(0)[i]*uNodes[i].getValue();
          Ad1d1+=Ndd.col(0)[i]*toEigen(cps[i]);
          ad2d2+=Ndd.col(1)[i]*uNodes[i].getValue();
          Ad2d2+=Ndd.col(1)[i]*toEigen(cps[i]);
          ad1d2+=Ndd.col(2)[i]*uNodes[i].getValue();
          Ad1d2+=Ndd.col(2)[i]*toEigen(cps[i]);
        }

        Eigen::Matrix<ScalarType,2,2> b;
        b<< ad1d1.dot(a3),ad1d2.dot(a3),
            ad1d2.dot(a3),ad2d2.dot(a3);

        Eigen::Matrix<double,2,2> B;
        B<< Ad1d1.dot(A3),Ad1d2.dot(A3),
            Ad1d2.dot(A3),Ad2d2.dot(A3);
        assert(B2.isApprox(B));
        Eigen::Matrix<ScalarType,2,2> a;
        a<<a1.squaredNorm(), a1.dot(a2),
            a2.dot(a1),a2.squaredNorm();
        const auto JcontraInv= calculateContravariantBaseVectors(J,A3);
        const Eigen::Matrix<double,2,3> Jloc = Dune::orthonormalizeMatrixColumns(J.transpose()).transpose();
        const Eigen::Matrix<double,2,2> id =J*JcontraInv.transpose();
        assert(Dune::FloatCmp::eq(std::abs(id(0,0)),1.0));
        assert(Dune::FloatCmp::eq(std::abs(id(1,1)),1.0));
        assert(Dune::FloatCmp::lt(std::abs(id(1,0)),1e-8));
        const auto epsV= toVoigt((0.5*(toLocalCartesian(a-A,JcontraInv,Jloc))).eval());
        const auto kappaV= toVoigt(toLocalCartesian(-b +B,JcontraInv,Jloc).eval());

        energy += 0.5*(epsV.dot(Cmat*epsV)*thickness_+kappaV.dot(Cmat*kappaV)*Dune::power(thickness_,3)/12.0) * geo.integrationElement(gp.position()) * gp.weight();
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

      const auto& element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary or not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine
            = Dune::QuadratureRules<double, Traits::mydim - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, Traits::mydim>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double intElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto u = uFunction.evaluate(quadPos);

          // Value of the Neumann data at the current position
          const auto neumannValue
              = neumannBoundaryLoad(toEigen(intersection.geometry().global(curQuad.position())), lambda);
          energy -= neumannValue.dot(u) * curQuad.weight() * intElement;
        }
      }

      return energy;
    }
  };
}  // namespace Ikarus
