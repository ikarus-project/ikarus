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
#include <ikarus/utils/observer/controlVTKWriter.hh>

namespace Ikarus {

  template <typename Basis_, typename Material_, typename FErequirements_ = FErequirements<>, bool useEigenRef = false>
  class NonLinearElastic : public PowerBasisFE<typename Basis_::FlatBasis>,
                           public Ikarus::AutoDiffFE<NonLinearElastic<Basis_, Material_, FErequirements_, useEigenRef>,
                                                     typename Basis_::FlatBasis, FErequirements_, useEigenRef> {
  public:
    using Basis             = Basis_;
    using FlatBasis         = typename Basis::FlatBasis;
    using BaseDisp          = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using BaseAD            = Ikarus::AutoDiffFE<NonLinearElastic<Basis_, Material_, FErequirements_, useEigenRef>,
                                      typename Basis_::FlatBasis, FErequirements_, useEigenRef>;
    using LocalView         = typename FlatBasis::LocalView;
    using Geometry          = typename LocalView::Element::Geometry;
    using GridView          = typename FlatBasis::GridView;
    using FERequirementType = FErequirements_;

    using ResultRequirementsType = ResultRequirements<FERequirementType>;

    using BaseAD::localView;
    using BaseAD::size;
    friend BaseAD;
    using Element  = typename LocalView::Element;
    using Material = Material_;

    template <typename VolumeLoad = LoadDefault, typename NeumannBoundaryLoad = LoadDefault>
    NonLinearElastic(const Basis& globalBasis, const typename LocalView::Element& element, const Material& p_mat,
                     VolumeLoad p_volumeLoad = {}, const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
                     NeumannBoundaryLoad p_neumannBoundaryLoad = {})
        : BaseDisp(globalBasis.flat(), element),
          BaseAD(globalBasis.flat(), element),
          neumannBoundary{p_neumannBoundary},
          mat{p_mat} {
      this->localView().bind(element);
      const int order = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      localBasis      = Dune::CachedLocalBasis(this->localView().tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(this->localView().element().type(), order),
                      Dune::bindDerivatives(0, 1));

      if constexpr (!std::is_same_v<VolumeLoad, LoadDefault>) volumeLoad = p_volumeLoad;
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, LoadDefault>) neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

    using Traits = TraitsFromLocalView<LocalView, useEigenRef>;

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
        const auto H        = uFunction.evaluateDerivative(gpIndex, Dune::wrt(spatialAll), Dune::on(gridElement));
        const auto E        = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
        const auto EVoigt   = toVoigt(E);
        auto internalEnergy = matAD.template storedEnergy<StrainTags::greenLagrangian>(EVoigt);
        energy += internalEnergy * geo.integrationElement(gp.position()) * gp.weight();
      }

      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
          const auto u                                 = uFunction.evaluate(gpIndex);
          Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigen(gp.position()), lambda);
          energy -= u.dot(fext) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      const int order = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary) return energy;

      auto element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary or not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, Traits::mydim - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, Traits::mydim>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto u = uFunction.evaluate(quadPos);

          // Value of the Neumann data at the current position
          auto neumannValue = neumannBoundaryLoad(toEigen(intersection.geometry().global(curQuad.position())), lambda);

          energy -= neumannValue.dot(u) * curQuad.weight() * integrationElement;
        }
      }

      return energy;
    }

  public:
    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      const auto& d      = req.getGlobalSolution(Ikarus::FESolutions::displacement);
      const auto& lambda = req.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      Dune::BlockVector<RealTuple<double, Traits::dimension>> disp(fe.size());

      for (auto i = 0U; i < fe.size(); ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      const auto geo = this->localView().element().geometry();

      Dune::StandardLocalFunction uFunction(localBasis, disp, std::make_shared<const Geometry>(geo));
      const auto H      = uFunction.evaluateDerivative(local, Dune::wrt(spatialAll), Dune::on(gridElement));
      const auto E      = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
      const auto EVoigt = toVoigt(E);
      auto PK2          = mat.template stresses<StrainTags::greenLagrangian>(EVoigt);

      typename ResultTypeMap<double>::ResultArray resultVector;
      if (req.isResultRequested(ResultType::PK2Stress)) {
        resultVector.resizeLike(PK2);
        resultVector = PK2;
        result.insertOrAssignResult(ResultType::PK2Stress, resultVector);
      }
    }

  private:
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
    Material mat;
  };

}  // namespace Ikarus
