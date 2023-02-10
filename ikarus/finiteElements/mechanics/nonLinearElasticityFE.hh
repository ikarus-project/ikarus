// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>
#include <iosfwd>

#include <dune/common/classname.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feBases/powerBasisFE.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename Basis>
  class NonLinearElasticityFE : public PowerBasisFE<Basis>,
                                public Ikarus::AutoDiffFE<NonLinearElasticityFE<Basis>, Basis> {
  public:
    using BaseDisp = PowerBasisFE<Basis>;  // Handles globalIndices function
    using BaseAD   = AutoDiffFE<NonLinearElasticityFE<Basis>, Basis>;
    using BaseAD::size;
    using GlobalIndex = typename PowerBasisFE<Basis>::GlobalIndex;
    friend BaseAD;
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    using LocalView         = typename Basis::LocalView;
    using Geometry          = typename LocalView::Element::Geometry;
    using GridView          = typename Basis::GridView;

    template <typename VolumeLoad, typename NeumannBoundaryLoad>
    NonLinearElasticityFE(Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu,
                          const BoundaryPatch<GridView>* neumannBoundary,
                          const NeumannBoundaryLoad& neumannBoundaryLoad, const VolumeLoad& p_volumeLoad)
        : BaseDisp(globalBasis, element),
          BaseAD(globalBasis, element),
          localView_{globalBasis.localView()},
          volumeLoad(p_volumeLoad),
          neumannBoundaryLoad_{neumannBoundaryLoad},
          neumannBoundary_{neumannBoundary},
          emod_{emod},
          nu_{nu} {
      localView_.bind(element);
      const int order = 2 * (localView_.tree().child(0).finiteElement().localBasis().order());
      localBasis      = Dune::CachedLocalBasis(localView_.tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),
                      Dune::bindDerivatives(0, 1));
    }

    using Traits = TraitsFromLocalView<LocalView>;

  private:
    template <class ScalarType>
    ScalarType calculateScalarImpl(const FERequirementType& req, Eigen::VectorX<ScalarType>& dx) const {
      const auto& d      = req.getGlobalSolution(Ikarus::FESolutions::displacement);
      const auto& lambda = req.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Dune::BlockVector<RealTuple<ScalarType, Traits::dimension>> disp(fe.size());

      for (auto i = 0U; i < fe.size(); ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp[i][k2] = dx[i * 2 + k2] + d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];

      ScalarType energy = 0.0;
      Eigen::Matrix3<ScalarType> C;
      C.setZero();  // plane stress
      C(0, 0) = C(1, 1) = 1;
      C(0, 1) = C(1, 0) = nu_;
      C(2, 2)           = (1 - nu_) / 2;
      C *= emod_ / (1 - nu_ * nu_);
      const auto geo = localView_.element().geometry();
      Dune::StandardLocalFunction uFunction(localBasis, disp, std::make_shared<const Geometry>(geo));
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto u      = uFunction.evaluate(gpIndex);
        const auto H      = uFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll), Dune::on(gridElement));
        const auto E      = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
        const auto EVoigt = toVoigt(E);

        Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigen(gp.position()), lambda);
        energy += (0.5 * EVoigt.dot(C * EVoigt) - u.dot(fext)) * geo.integrationElement(gp.position()) * gp.weight();
      }
      const int order = 2 * (localView_.tree().child(0).finiteElement().localBasis().order());
      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary_) return energy;

      auto element = localView_.element();
      for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
        if (not neumannBoundary_ or not neumannBoundary_->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, Traits::mydim - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, Traits::mydim>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto u = uFunction.evaluate(quadPos);

          // Value of the Neumann data at the current position
          auto neumannValue = neumannBoundaryLoad_(toEigen(intersection.geometry().global(curQuad.position())), lambda);

          energy -= neumannValue.dot(u) * curQuad.weight() * integrationElement;
        }
      }

      return energy;
    }

    LocalView localView_;
    Dune::CachedLocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        volumeLoad;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                          const double&)>
        neumannBoundaryLoad_;
    const BoundaryPatch<GridView>* neumannBoundary_;
    double emod_;
    double nu_;
  };

}  // namespace Ikarus
