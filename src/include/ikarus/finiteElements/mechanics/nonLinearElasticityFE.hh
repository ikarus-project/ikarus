
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
#include <ikarus/localFunctions/impl/standardLocalFunction.hh>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <concepts>
#include <iostream>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/interface/feTraits.hh>
#include <ikarus/finiteElements/interface/finiteElementFunctionConcepts.hh>
#include <ikarus/finiteElements/interface/interfaceFiniteElement.hh>
#include <ikarus/finiteElements/mechanics/displacementFE.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>

namespace Ikarus {

  template <typename Basis>
  class NonLinearElasticityFE
      : public DisplacementFE<Basis>,
        public Ikarus::AutoDiffFE<NonLinearElasticityFE<Basis>, Basis> {
  public:
    using BaseDisp = DisplacementFE<Basis>;  // Handles globalIndices function
    using BaseAD   = AutoDiffFE<NonLinearElasticityFE<Basis>, Basis>;
    using BaseAD::size;
    using GlobalIndex = typename DisplacementFE<Basis>::GlobalIndex;
    friend BaseAD;
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    using LocalView         = typename Basis::LocalView;

    template <typename VolumeLoad>
    NonLinearElasticityFE(Basis& globalBasis, const typename LocalView::Element& element, double emod,
                          double nu, const VolumeLoad& p_volumeLoad)
        : BaseDisp(globalBasis, element),
          BaseAD(globalBasis, element),
          localView_{globalBasis.localView()},
          volumeLoad(p_volumeLoad),
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
      Eigen::Matrix3<ScalarType> C;
      C.setZero();  // plane stress
      C(0, 0) = C(1, 1) = 1;
      C(0, 1) = C(1, 0) = nu_;
      C(2, 2)           = (1 - nu_) / 2;
      C *= emod_ / (1 - nu_ * nu_);
      const auto geo = localView_.element().geometry();
      Ikarus::StandardLocalFunction uFunction(localBasis, disp);
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto Jinv      = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().inverse().eval();
        const auto u      = uFunction.evaluateFunction(gpIndex);
        const auto H      = uFunction.evaluateDerivative(gpIndex, wrt(DerivativeDirections::spatialAll),
                                                    transformWith(Jinv));
        const auto E      = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
        const auto EVoigt = toVoigt(E);

        Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigenVector(gp.position()), lambda);
        energy += (0.5 * EVoigt.dot(C * EVoigt) - u.dot(fext)) * geo.integrationElement(gp.position()) * gp.weight();
      }
      return energy;
    }

    LocalView localView_;
    Ikarus::LocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        volumeLoad;
    double emod_;
    double nu_;
  };

}  // namespace Ikarus::FiniteElements
