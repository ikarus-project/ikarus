// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>
#include <iosfwd>
#include <optional>
#include <type_traits>

#include <dune/common/classname.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/localfefunctions/expressions/linearStrainsExpr.hh>
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
#include <ikarus/utils/defaultFunctions.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename Basis_, typename FErequirements_ = FErequirements<>, bool useEigenRef = false>
  class LinearElastic : public PowerBasisFE<typename Basis_::FlatBasis> {
  public:
    using Basis                  = Basis_;
    using FlatBasis              = typename Basis::FlatBasis;
    using BaseDisp               = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType      = FErequirements_;
    using ResultRequirementsType = ResultRequirements<FERequirementType>;
    using LocalView              = typename FlatBasis::LocalView;
    using Element                = typename LocalView::Element;
    using GridView               = typename FlatBasis::GridView;

    using Traits = TraitsFromLocalView<LocalView, useEigenRef>;

    static constexpr int mydim = Traits::mydim;

    template <typename VolumeLoad = LoadDefault, typename NeumannBoundaryLoad = LoadDefault>
    LinearElastic(const Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu,
                  VolumeLoad p_volumeLoad = {}, const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
                  NeumannBoundaryLoad p_neumannBoundaryLoad = {})
        : BaseDisp(globalBasis.flat(), element), neumannBoundary{p_neumannBoundary}, emod_{emod}, nu_{nu} {
      this->localView().bind(element);
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes     = fe.size();
      dispAtNodes.resize(fe.size());
      const int order = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      localBasis      = Dune::CachedLocalBasis(this->localView().tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(this->localView().element().type(), order),
                      Dune::bindDerivatives(0, 1));

      if constexpr (!std::is_same_v<VolumeLoad, LoadDefault>) volumeLoad = p_volumeLoad;
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, LoadDefault>) neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

  public:
    //    const auto& localView() const { return localView(); }

    auto getDisplacementFunction(const FERequirementType& par) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);

      for (auto i = 0U; i < dispAtNodes.size(); ++i)
        for (auto k2 = 0U; k2 < mydim; ++k2)
          dispAtNodes[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      auto geo = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          this->localView().element().geometry());
      Dune::StandardLocalFunction uFunction(localBasis, dispAtNodes, geo);

      return uFunction;
    }

    auto getStrainFunction(const FERequirementType& par) const { return linearStrains(getDisplacementFunction(par)); }

    auto getMaterialTangent() const {
      if constexpr (mydim == 2)
        return planeStressLinearElasticMaterialTangent(emod_, nu_);
      else if constexpr (mydim == 3)
        return linearElasticMaterialTangent3D(emod_, nu_);
    }

    auto getMaterialTangentFunction([[maybe_unused]] const FERequirementType& par) const {
      return [&]([[maybe_unused]] auto gp) { return getMaterialTangent(); };
    }

    double calculateScalar(const FERequirementType& par) const {
      const auto u       = getDisplacementFunction(par);
      const auto eps     = getStrainFunction(par);
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto C = getMaterialTangent();

      const auto geo = this->localView().element().geometry();
      double energy  = 0.0;
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const auto EVoigt = eps.evaluate(gpIndex, on(gridElement));

        energy += (0.5 * EVoigt.dot(C * EVoigt)) * geo.integrationElement(gp.position()) * gp.weight();
      }
      // External forces volume forces over the domain
      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
          const auto uVal                              = u.evaluate(gpIndex);
          Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigen(gp.position()), lambda);
          energy -= uVal.dot(fext) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary) return energy;

      auto element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, mydim - 1>::rule(intersection.type(), u.order());

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, mydim>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto uVal = u.evaluate(quadPos);

          // Value of the Neumann data at the current position
          auto neumannValue = neumannBoundaryLoad(toEigen(intersection.geometry().global(curQuad.position())), lambda);

          energy -= neumannValue.dot(uVal) * curQuad.weight() * integrationElement;
        }
      }
      return energy;
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType K) const {
      const auto eps = getStrainFunction(par);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto C   = getMaterialTangent();
      const auto geo = this->localView().element().geometry();

      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          for (size_t j = 0; j < numberOfNodes; ++j) {
            const auto bopJ = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), on(gridElement));
            K.template block<mydim, mydim>(i * mydim, j * mydim) += bopI.transpose() * C * bopJ * intElement;
          }
        }
      }
    }

    void calculateAt(const ResultRequirementsType& req, const Eigen::Vector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::Indices;
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto eps = getStrainFunction(req.getFERequirements());
      const auto C   = getMaterialTangent();
      auto gp        = toDune(local);
      auto epsVoigt  = eps.evaluate(gp, on(gridElement));

      auto cauchyStress = (C * epsVoigt).eval();

      typename ResultTypeMap<double>::ResultArray resultVector;
      if (req.isResultRequested(ResultType::linearStress)) {
        resultVector.resize(3, 1);
        resultVector = cauchyStress;
        result.insertOrAssignResult(ResultType::linearStress, resultVector);
      }
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType force) const {
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto eps     = getStrainFunction(par);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto C   = getMaterialTangent();
      const auto geo = this->localView().element().geometry();

      // Internal forces
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        auto stresses           = (C * eps.evaluate(gpIndex, on(gridElement))).eval();
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          force.template segment<mydim>(mydim * i) += bopI.transpose() * stresses * intElement;
        }
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        const auto u = getDisplacementFunction(par);
        for (const auto& [gpIndex, gp] : u.viewOverIntegrationPoints()) {
          Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigen(gp.position()), lambda);
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<mydim>(mydim * i)
                -= udCi * fext * geo.integrationElement(gp.position()) * gp.weight();
          }
        }
      }

      // External forces, boundary forces, i.e. at the Neumann boundary
      if (not neumannBoundary) return;

      const auto u = getDisplacementFunction(par);
      auto element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary->contains(intersection)) continue;

        // Integration rule along the boundary
        const auto& quadLine = Dune::QuadratureRules<double, mydim - 1>::rule(intersection.type(), u.order());

        for (const auto& curQuad : quadLine) {
          const Dune::FieldVector<double, mydim>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function wrt the i-th coef
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(quadPos, wrt(coeff(i)));

            // Value of the Neumann data at the current position
            auto neumannValue
                = neumannBoundaryLoad(toEigen(intersection.geometry().global(curQuad.position())), lambda);
            force.template segment<mydim>(mydim * i) -= udCi * neumannValue * curQuad.weight() * integrationElement;
          }
        }
      }
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
    double emod_;
    double nu_;
    size_t numberOfNodes{0};
  };

}  // namespace Ikarus
