
// /*
//  *  This file is part of the Ikarus distribution (https://github.com/rath3t/Ikarus).
//  *  Copyright (c) 2021 Alexander Müller.
//  *  Institut fuer Baustatik und Baudynamik
//  *  Universität Stuttgart
//  *
//  *  This library is free software; you can redistribute it and/or
//  *   modify it under the terms of the GNU Lesser General Public
//  *   License as published by the Free Software Foundation; either
//  *   version 2.1 of the License, or (at your option) any later version.
//
// *   This library is distributed in the hope that it will be useful,
// *   but WITHOUT ANY WARRANTY; without even the implied warranty of
// *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// *   Lesser General Public License for more details.
//
// *   You should have received a copy of the GNU Lesser General Public
// *   License along with this library; if not, write to the Free Software
// *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// *  USA
// *

#pragma once
#include "src/include/ikarus/finiteElements/feTraits.hh"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <concepts>
#include <iosfwd>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/fufem/boundarypatch.hh>

#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/mechanics/displacementFE.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/localFunctions/impl/standardLocalFunction.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename Basis>
  class Q1LinearElastic : public DisplacementFE<Basis>{
  public:
    using BaseDisp = DisplacementFE<Basis>;  // Handles globalIndices function
//    using BaseAD   = AutoDiffFE<Q1LinearElastic<Basis>, Basis>;
//    using BaseAD::size;
    using GlobalIndex = typename DisplacementFE<Basis>::GlobalIndex;
//    friend BaseAD;
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    using LocalView         = typename Basis::LocalView;
    using GridView         = typename Basis::GridView;

    template <typename VolumeLoad, typename NeumannBoundaryLoad>
    Q1LinearElastic(Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu, const BoundaryPatch<GridView> * neumannBoundary,
                    const NeumannBoundaryLoad& neumannBoundaryLoad,
                    const VolumeLoad& p_volumeLoad)
        : BaseDisp(globalBasis, element),
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

  public:
    double calculateScalar(const FERequirementType& par)const {
      const auto& d      = par.getSolution(Ikarus::FESolutions::displacement);
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Dune::BlockVector<Ikarus::RealTuple<double, Traits::dimension>> disp(fe.size());

      for (auto i = 0U; i < fe.size(); ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp[i][k2] = d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];

      double energy = 0.0;
      const int order   = 2 * (fe.localBasis().order());
      const auto& rule  = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
      Eigen::Matrix3<double> C = planeStressLinearElasticMaterialTangent(emod_,nu_);
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

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      const auto& d      = par.getSolution(Ikarus::FESolutions::displacement);
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace DerivativeDirections;

      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Dune::BlockVector<Ikarus::RealTuple<double, Traits::dimension>> disp(fe.size());

      for (auto i = 0U; i < fe.size(); ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp[i][k2] = d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];

      const int order   = 2 * (fe.localBasis().order());
      const auto& rule  = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
      Eigen::Matrix3<double> C = planeStressLinearElasticMaterialTangent(emod_,nu_);
      const auto geo = localView_.element().geometry();
      Ikarus::StandardLocalFunction uFunction(localBasis, disp);
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto Jinv = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().inverse().eval();
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        for (size_t i = 0; i < fe.size(); ++i) {
          const auto dHdCi
              = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i)), transformWith(Jinv));
          Eigen::Matrix<double,3,2> bopI;
          bopI.row(0)<< dHdCi[0].diagonal()(0),0;
          bopI.row(1)<< 0,dHdCi[1].diagonal()(1);
          bopI.row(2)<< dHdCi[1].diagonal()(0),dHdCi[0].diagonal()(1);
          for (size_t j = 0; j < fe.size(); ++j) {
            const auto dHdCj
                = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(j)), transformWith(Jinv));
            Eigen::Matrix<double,3,2> bopJ;
            bopJ.row(0)<< dHdCj[0].diagonal()(0),0;
            bopJ.row(1)<< 0,dHdCj[1].diagonal()(1);
            bopJ.row(2)<< dHdCj[1].diagonal()(0),dHdCj[0].diagonal()(1);
            h.template block<2, 2>(i*Traits::mydim,j*Traits::mydim)+= bopI.transpose()*C*bopJ*intElement;
          }
        }
      }
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      const auto& d      = par.getSolution(Ikarus::FESolutions::displacement);
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace DerivativeDirections;

      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Dune::BlockVector<Ikarus::RealTuple<double, Traits::dimension>> disp(fe.size());

      for (auto i = 0U; i < fe.size(); ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp[i][k2] = d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];

      const int order   = 2 * (fe.localBasis().order());
      const auto& rule  = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
      Eigen::Matrix3<double> C = planeStressLinearElasticMaterialTangent(emod_,nu_);
      const auto geo = localView_.element().geometry();
      Ikarus::StandardLocalFunction uFunction(localBasis, disp);
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto Jinv = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().inverse().eval();
        const auto u    = uFunction.evaluateFunction(gpIndex);
        Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigenVector(gp.position()), lambda);
        for (size_t i = 0; i < fe.size(); ++i) {
//          const auto gradUdCi =uFunction.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i)), transformWith(Jinv));
          const auto udCi =uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));

          g.template segment<Traits::mydim>(Traits::mydim*i) -= udCi*fext* geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      //line or surface loads, i.e. neumann boundary
      if (not neumannBoundary_) return;

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

          // The value of the local function wrt the i-th coef
          for (size_t i = 0; i < fe.size(); ++i) {
            const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

            // Value of the Neumann data at the current position
            auto neumannValue
                = neumannBoundaryLoad_(toEigenVector(intersection.geometry().global(curQuad.position())), lambda);
            g.template segment<Traits::mydim>(Traits::mydim * i)
                -= udCi * neumannValue * curQuad.weight() * integrationElement;
          }
        }
      }
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

}  // namespace Ikarus
