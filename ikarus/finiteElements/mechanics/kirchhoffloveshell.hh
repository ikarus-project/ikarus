// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "fesettings.hh"

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
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/utils/tensorUtils.hh>
namespace Ikarus {

  template <typename Basis_, typename FERequirements_ = FErequirements<>, bool useEigenRef = false>
  class KirchhoffLoveShell : public PowerBasisFE<typename Basis_::FlatBasis> {
  public:
    using Basis                   = Basis_;
    using FlatBasis               = typename Basis::FlatBasis;
    using BasePowerFE             = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType       = FERequirements_;
    using ResultRequirementsType  = ResultRequirements<FERequirementType>;
    using LocalView               = typename FlatBasis::LocalView;
    using Element                 = typename LocalView::Element;
    using Geometry                = typename Element::Geometry;
    using GridView                = typename FlatBasis::GridView;
    using Traits                  = TraitsFromLocalView<LocalView, useEigenRef>;
    static constexpr int myDim    = Traits::mydim;
    static constexpr int worlddim = Traits::worlddim;

    template <typename VolumeLoad = LoadDefault, typename NeumannBoundaryLoad = LoadDefault>
    KirchhoffLoveShell(const Basis& globalBasis, const typename LocalView::Element& element,
                       const FESettings& p_feSettings, VolumeLoad p_volumeLoad = {},
                       const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
                       NeumannBoundaryLoad p_neumannBoundaryLoad        = {})
        : BasePowerFE(globalBasis.flat(), element), neumannBoundary{p_neumannBoundary}, fESettings{p_feSettings} {
      this->localView().bind(element);
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes     = fe.size();
      dispAtNodes.resize(fe.size());
      dispAtNodesDual.resize(fe.size());
      dispAtNodesDual2nd.resize(fe.size());
      order      = 2 * (fe.localBasis().order());
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
    template <typename ScalarType = double>
    auto getDisplacementFunction(const FERequirementType& par,
                                 const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);
      auto geo      = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          this->localView().element().geometry());
      //      Dune::BlockVector<Dune::RealTuple<ScalarType, Traits::worlddim>> disp(dispAtNodes.size());
      if constexpr (std::is_same_v<ScalarType, autodiff::dual>) {
        if (dx)
          for (auto i = 0U; i < dispAtNodesDual.size(); ++i)
            for (auto k2 = 0U; k2 < worlddim; ++k2)
              dispAtNodesDual[i][k2]
                  = dx.value()[i * worlddim + k2]
                    + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
        else
          for (auto i = 0U; i < dispAtNodesDual.size(); ++i)
            for (auto k2 = 0U; k2 < worlddim; ++k2)
              dispAtNodesDual[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      } else if constexpr (std::is_same_v<ScalarType, autodiff::dual2nd>) {
        if (dx)
          for (auto i = 0U; i < dispAtNodesDual2nd.size(); ++i)
            for (auto k2 = 0U; k2 < worlddim; ++k2)
              dispAtNodesDual2nd[i][k2]
                  = dx.value()[i * worlddim + k2]
                    + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
        else
          for (auto i = 0U; i < dispAtNodesDual2nd.size(); ++i)
            for (auto k2 = 0U; k2 < worlddim; ++k2)
              dispAtNodesDual2nd[i][k2]
                  = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      } else if constexpr (std::is_same_v<ScalarType, double>) {
        if (dx)
          for (auto i = 0U; i < dispAtNodes.size(); ++i)
            for (auto k2 = 0U; k2 < worlddim; ++k2)
              dispAtNodes[i][k2] = dx.value()[i * worlddim + k2]
                                   + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
        else
          for (auto i = 0U; i < dispAtNodes.size(); ++i)
            for (auto k2 = 0U; k2 < worlddim; ++k2)
              dispAtNodes[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      }

      if constexpr (std::is_same_v<ScalarType, autodiff::dual>) {
        Dune::StandardLocalFunction uFunction(localBasis, dispAtNodesDual, geo);
        return std::make_pair(uFunction, std::ref(dispAtNodesDual));
      } else if constexpr (std::is_same_v<ScalarType, autodiff::dual2nd>) {
        Dune::StandardLocalFunction uFunction(localBasis, dispAtNodesDual2nd, geo);
        return std::make_pair(uFunction, std::ref(dispAtNodesDual2nd));
      } else if constexpr (std::is_same_v<ScalarType, double>) {
        Dune::StandardLocalFunction uFunction(localBasis, dispAtNodes, geo);
        return std::make_pair(uFunction, std::ref(dispAtNodes));
      }
    }

    inline double calculateScalar(const FERequirementType& par) const { return calculateScalarImpl<double>(par); }

    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      DUNE_THROW(Dune::NotImplemented, "No results are implemented");
    }

    inline void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(par, force);
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
    mutable Dune::BlockVector<Dune::RealTuple<double, Traits::worlddim>> dispAtNodes;
    mutable Dune::BlockVector<Dune::RealTuple<autodiff::dual, Traits::worlddim>> dispAtNodesDual;
    mutable Dune::BlockVector<Dune::RealTuple<autodiff::dual2nd, Traits::worlddim>> dispAtNodesDual2nd;
    FESettings fESettings;
    size_t numberOfNodes{0};
    int order{};

  protected:
    mutable std::vector<Dune::FieldVector<double, 2>> lagrangePoints;
    mutable std::vector<double> out;
    template <typename ScalarType>
    auto calculateScalarImpl(const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                                           = std::nullopt) const -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [uFunction, uNodes] = getDisplacementFunction(par, dx);
      const auto& lambda             = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto geo                 = this->localView().element().geometry();
      ScalarType energy              = 0.0;
      const auto uasMatrix           = Dune::viewAsEigenMatrixAsDynFixed(uNodes);

      const auto& simFlag = fESettings.request<int>("simulationFlag");

      Dune::LagrangeCubeLocalFiniteElement<double, double, 2, 2> q1lfem2D;

      if (simFlag == 1) {
        lagrangePoints.resize(q1lfem2D.size());
        for (int i = 0; i < 2; i++) {
          auto ithCoord = [&i](const FieldVector<double, 2>& x) { return x[i]; };
          q1lfem2D.localInterpolation().interpolate(ithCoord, out);
          for (std::size_t jI = 0; jI < out.size(); jI++)
            lagrangePoints[jI][i] = out[jI];
        }
      } else if (simFlag == 2) {
        lagrangePoints.resize(q1lfem2D.size());
        lagrangePoints[0] = {0, 0.5};
        lagrangePoints[1] = {0.5, 0};
        lagrangePoints[2] = {1, 0.5};
        lagrangePoints[3] = {0.5, 1};
      }

      std::vector<Eigen::Vector<ScalarType, 3>> membraneStrainsAtVertices;
      for (int i = 0; auto& lP : lagrangePoints) {
        const auto J                                = toEigen(geo.jacobianTransposed(lP));
        const Eigen::Matrix<double, 2, 2> A         = J * J.transpose();
        const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
            uFunction.evaluateDerivative(lP, wrt(spatialAll, Dune::on(DerivativeDirections::referenceElement))));
        // const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();

        const auto  epsV = toVoigt((0.5 * (gradu.transpose()*J.transpose() + J*gradu +gradu.transpose()*gradu)).eval()).eval();
        membraneStrainsAtVertices.push_back(epsV);
      }

      std::vector<Dune::FieldVector<double, 1>> NANS;

      auto evalMembraneStrains = [&](auto& gppos) {
        if (simFlag == 1)
          q1lfem2D.localBasis().evaluateFunction(gppos, NANS);
        else if (simFlag == 2) {
          NANS.resize(4);
          NANS[0] = 0.5 * (1 - gppos[1]);
          NANS[1] = 0.5 * (1 - gppos[0]);
          NANS[2] = 0.5 * (1 + gppos[1]);
          NANS[3] = 0.5 * (1 + gppos[0]);
        }
        Eigen::Vector<ScalarType, 3> res;
        res.setZero();
        for (int i = 0; i < NANS.size(); ++i) {
          res += membraneStrainsAtVertices[i] * NANS[i][0];
        }
        return res;
      };

      //      auto evalMembraneStrains = [&](auto& gppos)
      //      {
      //        std::vector<Dune::FieldVector<double,1>> N;
      //        q1lfem2D.localBasis().evaluateFunction(gppos,N);
      //        Eigen::Vector<ScalarType,3> res;
      //        res.setZero();
      //        for (int i = 0; i < 4; ++i) {
      //          res+=membraneStrainsAtVertices[i]*N[i][0];
      //        }
      //        return res;
      //      };

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto [X, Jd, Hd]                      = geo.impl().zeroFirstAndSecondDerivativeOfPosition(gp.position());
        const auto J                                = toEigen(Jd);
        const auto H                                = toEigen(Hd);
        const Eigen::Matrix<double, 2, 2> A         = J * J.transpose();
        const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
            uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, Dune::on(DerivativeDirections::referenceElement))));
        const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();

        const auto& Ndd                     = localBasis.evaluateSecondDerivatives(gpIndex);
        const auto h                        = (H + Ndd.transpose().template cast<ScalarType>() * uasMatrix).eval();
        const Eigen::Vector3<ScalarType> a3 = (j.row(0).cross(j.row(1))).normalized();
        Eigen::Vector<ScalarType, 3> bV     = h * a3;
        bV(2) *= 2;  // Voigt notation requires the two here

        Eigen::Matrix<double, 3, 3> G;
        G.setZero();
        G.block<2, 2>(0, 0)                    = A;
        G(2, 2)                                = 1;
        const Eigen::Matrix<double, 3, 3> GInv = G.inverse();

        const auto& thickness_ = fESettings.request<double>("thickness");

        Eigen::Vector3<ScalarType> epsV;
        if (simFlag == 0)          //2x3  2x3                    2x3 3x2
          epsV = toVoigt((0.5 * (gradu.transpose()*J.transpose() + J*gradu +gradu.transpose()*gradu)).eval()).eval();
        else if (simFlag == 1 or simFlag == 2)
          epsV = evalMembraneStrains(gp.position());

        //        const auto epsV = evalMembraneStrains(gp.position());

        // const auto epsV                 = toVoigt((0.5 * (j * j.transpose() - A)).eval()).eval();
        const auto BV                   = toVoigt(toEigen(geo.impl().secondFundamentalForm(gp.position())));
        const auto kappaV               = (BV - bV).eval();
        const auto C                    = materialTangent(GInv);
        const ScalarType membraneEnergy = 0.5 * thickness_ * epsV.dot(C * epsV);
        const ScalarType bendingEnergy  = 0.5 * Dune::power(thickness_, 3) / 12.0 * kappaV.dot(C * kappaV);
        energy += (membraneEnergy + bendingEnergy) * geo.integrationElement(gp.position()) * gp.weight();
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

        const auto& quadLine = Dune::QuadratureRules<double, Traits::mydim - 1>::rule(intersection.type(), order);

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

        template <typename ScalarType>
        void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                                 const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
          using namespace Dune::DerivativeDirections;
          using namespace Dune;
          const auto [uFunction, uNodes] = getDisplacementFunction(par, dx);
          const auto& lambda             = par.getParameter(Ikarus::FEParameter::loadfactor);
          const auto geo                 = this->localView().element().geometry();
          ScalarType energy              = 0.0;
          const auto uasMatrix           = Dune::viewAsEigenMatrixAsDynFixed(uNodes);
          const auto& thickness_ = fESettings.request<double>("thickness");

          // Internal forces
          for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
            const double intElement = geo.integrationElement(gp.position()) * gp.weight();
            const auto [X, Jd, Hd]                      =
            geo.impl().zeroFirstAndSecondDerivativeOfPosition(gp.position()); const auto J = toEigen(Jd); const auto H
            = toEigen(Hd); const Eigen::Matrix<double, 2, 2> A         = J * J.transpose(); const
            Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
                uFunction.evaluateDerivative(gpIndex, wrt(spatialAll,
                Dune::on(DerivativeDirections::referenceElement))));
            const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();

            const auto& Ndd                     = localBasis.evaluateSecondDerivatives(gpIndex);
            const auto& Nd                     = localBasis.evaluateJacobian(gpIndex);
            const auto h                        = H + Ndd.transpose().template cast<ScalarType>() * uasMatrix;
            const Eigen::Vector3<ScalarType> a3N = j.row(0).cross(j.row(1));
            const Eigen::Vector3<ScalarType> a3 = a3N.normalized();
            Eigen::Vector<ScalarType, 3> bV     = h * a3;
            bV(2) *= 2;  // Voigt notation requires the two here

            Eigen::Matrix<double, 3, 3> G;
            G.setZero();
            G.block<2, 2>(0, 0)                    = A;
            G(2, 2)                                = 1;
            const Eigen::Matrix<double, 3, 3> GInv = G.inverse();

//            const auto a3dI = a3

            const auto epsV                 = toVoigt((0.5 * (j * j.transpose() - A)).eval()).eval();
            const auto BV                   = toVoigt(toEigen(geo.impl().secondFundamentalForm(gp.position())));
            const auto kappaV               = (BV - bV).eval();
            const auto C = materialTangent(GInv);
            const Eigen::Vector<ScalarType,3> membraneForces = thickness_*C*epsV;
            const Eigen::Vector<ScalarType,3> moments = Dune::power(thickness_, 3) / 12.0*C*kappaV;

            for (size_t i = 0; i < numberOfNodes; ++i) {
              Eigen::Matrix<ScalarType, 3, 3> bopIMembrane = bopMembrane(j,Nd,i);
              force.template segment<worlddim>(worlddim * i) += bopIMembrane.transpose() * membraneForces * intElement;
            }
          }

          // External forces volume forces over the domain
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

          const auto u = getDisplacementFunction(par, dx);
          auto element = this->localView().element();
          for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
            if (not neumannBoundary->contains(intersection)) continue;

            // Integration rule along the boundary
            const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), order);

            for (const auto& curQuad : quadLine) {
              const Dune::FieldVector<double, myDim>& quadPos =
              intersection.geometryInInside().global(curQuad.position());

              const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

              // The value of the local function wrt the i-th coef
              for (size_t i = 0; i < numberOfNodes; ++i) {
                const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

                // Value of the Neumann data at the current position
                auto neumannValue
                    = neumannBoundaryLoad(toEigen(intersection.geometry().global(curQuad.position())), lambda);
                force.template segment<worlddim>(worlddim * i) -= udCi * neumannValue * curQuad.weight() *
                integrationElement;
              }
            }
          }
        }

  private:
    template <typename ScalarType>
    Eigen::Matrix<ScalarType, 3, 3> bopMembrane(const Eigen::Matrix<ScalarType, 2, 3>& J, const auto& dN,
                                                const int node)const  {
      Eigen::Matrix<ScalarType, 3, 3> bop;
      bop << J.row(0) * dN(node, 0), J.row(1) * dN(node, 1), J.row(0) * dN(node, 1) + J.row(1) * dN(node, 0);

      return bop;
    }

    Eigen::Matrix<double, 3, 3> materialTangent(const Eigen::Matrix<double, 3, 3>& Aconv) const {
      const auto& emod_     = fESettings.request<double>("youngs_modulus");
      const auto& nu_       = fESettings.request<double>("poissons_ratio");
      const double lambda   = emod_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
      const double mu       = emod_ / (2.0 * (1.0 + nu_));
      const double lambdbar = 2.0 * lambda * mu / (lambda + 2.0 * mu);
      Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>> moduli;
      const auto AconvT = TensorCast(Aconv, std::array<Eigen::Index, 2>({3, 3}));
      moduli = lambdbar * dyadic(AconvT, AconvT).eval() + 2 * mu * symmetricFourthOrder<double>(Aconv, Aconv);

      auto C   = Ikarus::toVoigt(moduli);
      Eigen::Matrix<double, 3, 3> C33 = C({0, 1, 5}, {0, 1, 5});
      return C33;
    }
  };
}  // namespace Ikarus
