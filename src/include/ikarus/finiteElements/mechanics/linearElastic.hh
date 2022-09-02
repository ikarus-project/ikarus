/*
* This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
* Copyright (c) 2022. The Ikarus developers.
*
* This library is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/

#pragma once
#include "src/include/ikarus/finiteElements/feTraits.hh"

#include <concepts>
#include <iosfwd>

#include <dune/common/classname.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feBases/powerBasisFE.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/localFunctions/expressions/linearStrainsExpr.hh>
#include <ikarus/localFunctions/impl/standardLocalFunction.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename Basis>
  class LinearElastic : public PowerBasisFE<Basis> {
  public:
    using BaseDisp          = PowerBasisFE<Basis>;  // Handles globalIndices function
    using GlobalIndex       = typename PowerBasisFE<Basis>::GlobalIndex;
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    using LocalView         = typename Basis::LocalView;
    using GridView          = typename Basis::GridView;

    template <typename VolumeLoad, typename NeumannBoundaryLoad>
    LinearElastic(Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu,
                  const BoundaryPatch<GridView>* neumannBoundary, const NeumannBoundaryLoad& neumannBoundaryLoad,
                  const VolumeLoad& p_volumeLoad)
        : BaseDisp(globalBasis, element),
          localView_{globalBasis.localView()},
          volumeLoad(p_volumeLoad),
          neumannBoundaryLoad_{neumannBoundaryLoad},
          neumannBoundary_{neumannBoundary},
          emod_{emod},
          nu_{nu} {
      localView_.bind(element);
      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes=fe.size();
      dispAtNodes.resize(fe.size());
      const int order = 2 * (localView_.tree().child(0).finiteElement().localBasis().order());
      localBasis      = Ikarus::LocalBasis(localView_.tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),
                      bindDerivatives(0, 1));
    }

    using Traits = TraitsFromLocalView<LocalView>;

    static constexpr int mydim = Traits::mydim;
    static constexpr int strainSize = mydim == 1 ? 1 : (mydim == 2 ? 3 : 6);

  public:
    const auto& getLocalView() const { return localView_; }

    auto getDisplacementFunction(const FERequirementType& par) const {
      const auto& d = par.getSolution(Ikarus::FESolutions::displacement);

      for (auto i = 0U; i < dispAtNodes.size(); ++i)
        for (auto k2 = 0U; k2 < mydim; ++k2)
          dispAtNodes[i][k2] = d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];

      Ikarus::StandardLocalFunction uFunction(localBasis, dispAtNodes);

      return uFunction;
    }

    auto getStrainFunction(const FERequirementType& par) const {
      return linearStrains(getDisplacementFunction(par));
    }

    auto getMaterialTangent() const {
        if constexpr (mydim == 2)
          return planeStressLinearElasticMaterialTangent(emod_, nu_);
        else if constexpr (mydim == 3)
          return LinearElasticMaterialTangent3D(emod_, nu_);
    }

    auto getMaterialTangentFunction(const FERequirementType& par) const {
      return [&](auto gp) {
        return getMaterialTangent();
      };
    }

    double calculateScalar(const FERequirementType& par) const {
      const auto u = getDisplacementFunction(par);
      const auto eps = getStrainFunction(par);
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

      const auto C = getMaterialTangent();

      const auto geo = localView_.element().geometry();
      double energy    = 0.0;
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const auto Jinv = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().inverse().eval();
        const auto EVoigt = eps.evaluateFunction(gpIndex, transformWith(Jinv));

        const auto uVal = u.evaluateFunction(gpIndex);
        Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigenVector(gp.position()), lambda);
        energy += (0.5 * EVoigt.dot(C * EVoigt) - uVal.dot(fext)) * geo.integrationElement(gp.position()) * gp.weight();
      }

      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary_) return energy;

      auto element = localView_.element();
      for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
        if (not neumannBoundary_ or not neumannBoundary_->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, mydim - 1>::rule(intersection.type(), u.order());

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, mydim>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto uVal = u.evaluateFunction(quadPos);

          // Value of the Neumann data at the current position
          auto neumannValue
              = neumannBoundaryLoad_(toEigenVector(intersection.geometry().global(curQuad.position())), lambda);

          energy -= neumannValue.dot(uVal) * curQuad.weight() * integrationElement;
        }
      }

      return energy;
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      const auto eps = getStrainFunction(par);
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace DerivativeDirections;

      const auto C = getMaterialTangent();
      const auto geo = localView_.element().geometry();

      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const auto Jinv         = toEigenMatrix(geo.jacobianTransposed(gp.position())).transpose().inverse().eval();
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), transformWith(Jinv));
          for (size_t j = 0; j < numberOfNodes; ++j) {
            const auto bopJ = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), transformWith(Jinv));
            h.template block<mydim, mydim>(i * mydim, j * mydim)
                += bopI.transpose() * C * bopJ * intElement;
          }
        }
      }
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace DerivativeDirections;

      const auto geo   = localView_.element().geometry();
       const auto u = getDisplacementFunction(par);
      for (const auto& [gpIndex, gp] : u.viewOverIntegrationPoints()) {
        Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigenVector(gp.position()), lambda);
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto udCi = u.evaluateDerivative(gpIndex, wrt(coeff(i)));
          g.template segment<mydim>(mydim * i)
              -= udCi * fext * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary_) return;

      auto element = localView_.element();
      for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
        if (not neumannBoundary_ or not neumannBoundary_->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, mydim - 1>::rule(intersection.type(), u.order());

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, mydim>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function wrt the i-th coef
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(quadPos, wrt(coeff(i)));

            // Value of the Neumann data at the current position
            auto neumannValue
                = neumannBoundaryLoad_(toEigenVector(intersection.geometry().global(curQuad.position())), lambda);
            g.template segment<mydim>(mydim * i)
                -= udCi * neumannValue * curQuad.weight() * integrationElement;
          }
        }
      }
    }

    LocalView localView_;
    Ikarus::LocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    // TODO: write as optional
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        volumeLoad;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        neumannBoundaryLoad_;
    const BoundaryPatch<GridView>* neumannBoundary_;
    mutable Dune::BlockVector<Ikarus::RealTuple<double, Traits::dimension>> dispAtNodes;
    double emod_;
    double nu_;
    size_t numberOfNodes{0};
  };

}  // namespace Ikarus
