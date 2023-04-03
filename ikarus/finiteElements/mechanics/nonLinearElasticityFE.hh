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
#include <ikarus/finiteElements/mechanics/materials.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename Basis, typename Material>
  class NonLinearElasticityFE
      : public PowerBasisFE<typename Basis::FlatBasis>,
        public Ikarus::AutoDiffFE<NonLinearElasticityFE<Basis, Material>, typename Basis::FlatBasis> {
  public:
    using FlatBasis = typename Basis::FlatBasis;
    using BaseDisp  = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using BaseAD    = Ikarus::AutoDiffFE<NonLinearElasticityFE<Basis, Material>, typename Basis::FlatBasis>;
    using BaseAD::localView;
    using BaseAD::size;
    friend BaseAD;
    using FERequirementType = FErequirements<Eigen::VectorXd>;
    using LocalView         = typename FlatBasis::LocalView;
    using Geometry          = typename LocalView::Element::Geometry;
    using GridView          = typename FlatBasis::GridView;

    template <typename VolumeLoad, typename NeumannBoundaryLoad>
    NonLinearElasticityFE(Basis& globalBasis, const typename LocalView::Element& element, const Material& p_mat,
                          const BoundaryPatch<GridView>* neumannBoundary,
                          const NeumannBoundaryLoad& neumannBoundaryLoad, const VolumeLoad& p_volumeLoad)
        : BaseDisp(globalBasis.flat(), element),
          BaseAD(globalBasis.flat(), element),
          volumeLoad(p_volumeLoad),
          neumannBoundaryLoad_{neumannBoundaryLoad},
          neumannBoundary_{neumannBoundary},
          mat{p_mat} {
      this->localView().bind(element);
      const int order = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      localBasis      = Dune::CachedLocalBasis(this->localView().tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(this->localView().element().type(), order),
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
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Dune::BlockVector<RealTuple<ScalarType, Traits::dimension>> disp(fe.size());

      for (auto i = 0U; i < fe.size(); ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp[i][k2]
              = dx[i * 2 + k2] + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];

      ScalarType energy = 0.0;

      decltype(auto) matAD = mat.template rebind<ScalarType>();

      const auto geo = this->localView().element().geometry();
      Dune::StandardLocalFunction uFunction(localBasis, disp, std::make_shared<const Geometry>(geo));
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto u        = uFunction.evaluate(gpIndex);
        const auto H        = uFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll), Dune::on(gridElement));
        const auto E        = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
        const auto EVoigt   = toVoigt(E);
        auto internalEnergy = matAD.template storedEnergy<StrainTags::greenLagrangian>(EVoigt);
        Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigen(gp.position()), lambda);
        energy += (internalEnergy - u.dot(fext)) * geo.integrationElement(gp.position()) * gp.weight();
      }
      const int order = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary_) return energy;

      auto element = this->localView().element();
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
    Material mat;
  };

}  // namespace Ikarus
