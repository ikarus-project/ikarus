
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
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <concepts>
#include <iostream>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include "ikarus/LocalBasis/localBasis.h"
#include "ikarus/utils/LinearAlgebraHelper.h"
#include <ikarus/FiniteElements/FEPolicies.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/FiniteElements/physicsHelper.h>
#include <ikarus/Geometries/GeometryWithExternalInput.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>
#include <ikarus/FiniteElements/AutodiffFE.h>

namespace Ikarus::FiniteElements {

  template <typename Basis>
  class NonLinearElasticityFEWithLocalBasis : public FEDisplacement<Basis>, public Ikarus::AutoDiffFEClean<NonLinearElasticityFEWithLocalBasis<Basis>, Basis> {
  public:
    using BaseDisp          = Ikarus::FiniteElements::FEDisplacement<Basis>; //Handles globalIndices function
    using BaseAD            = Ikarus::AutoDiffFEClean<NonLinearElasticityFEWithLocalBasis<Basis>, Basis>;
    friend BaseAD ;
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    using LocalView = typename Basis::LocalView;
    NonLinearElasticityFEWithLocalBasis(Basis& globalBasis,const typename LocalView::Element& element, double emod, double nu)
        :  BaseDisp(globalBasis,element), BaseAD(globalBasis,element), localView_{globalBasis.localView()}, emod_{emod}, nu_{nu} {

      localView_.bind(element);
      const int order  = 2 * (localView_.tree().child(0).finiteElement().localBasis().order());
      localBasis = Ikarus::LocalBasis(localView_.tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),0,1);
    }

    using Traits = TraitsFromLocalView<LocalView>;
    template <typename ST>
    using DefoGeo = Ikarus::Geometry::GeometryWithExternalInput<ST, Traits::mydim, Traits::dimension>;

  private:
    template <class ScalarType>
    ScalarType calculateScalarImpl(const FERequirementType& par,
                                   Eigen::VectorX<ScalarType>& dx) const {
      const auto& d      = par.sols[0].get();
      const auto& lambda = par.parameter.at(FEParameter::loadfactor);
      Eigen::VectorX<ScalarType> localDisp(localView_.size());
      localDisp.setZero();
      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Eigen::Matrix<ScalarType, Traits::dimension, Eigen::Dynamic> disp;
      disp.setZero(Eigen::NoChange, fe.size());
      for (auto i = 0U; i < fe.size(); ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp.col(i)(k2)
              = dx[i * 2 + k2] + d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];
      ScalarType energy = 0.0;

      const int order  = 2 * (fe.localBasis().order());
      const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
      Eigen::Matrix3<ScalarType> C;
      C.setZero();
      C(0, 0) = C(1, 1) = 1;
      C(0, 1) = C(1, 0) = nu_;
      C(2, 2)           = (1 - nu_) / 2;
      C *= emod_ / (1 - nu_ * nu_);
      const auto geo = localView_.element().geometry();
//      Ikarus::LocalBasis localBasis(fe.localBasis());
//      Eigen::Matrix<double, Eigen::Dynamic, Traits::mydim> dN;
//      Eigen::VectorXd N;
      for (const auto& [index, gp ,N, dN] : localBasis.viewOverFunctionAndJacobian()) {
        const auto J = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().eval();
//        localBasis.evaluateFunctionAndJacobian(gp.position(), N, dN);
        const Eigen::Vector<double, Traits::worlddim> X = toEigenVector(geo.global(gp.position()));
        Eigen::Vector<ScalarType, Traits::worlddim> x   = X;
        for (int i = 0; i < N.size(); ++i)
          x += disp.col(i) * N[i];

        const auto H      = DefoGeo<ScalarType>::jacobianTransposed(dN*J.inverse(), disp).eval();
        const auto E      = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
        const auto EVoigt = toVoigt(E);

        Eigen::Vector<double, Traits::worlddim> fext;
        fext.setZero();
        fext[1] = 2*lambda;
        fext[0] = lambda;
        energy += (0.5*EVoigt.dot(C * EVoigt) - x.dot(fext)) * geo.integrationElement(gp.position()) * gp.weight();
      }
      return energy;
    }

    LocalView localView_;
    Ikarus::LocalBasis< std::remove_cvref_t<decltype( std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>> localBasis;
    double emod_;
    double nu_;
  };

}  // namespace Ikarus::FiniteElements
