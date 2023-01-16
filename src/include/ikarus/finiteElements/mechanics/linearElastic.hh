// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>
#include <iosfwd>
#include <optional>

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
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename Basis>
  class LinearElastic : public PowerBasisFE<Basis> {
  public:
    using BaseDisp               = PowerBasisFE<Basis>;  // Handles globalIndices function
    using GlobalIndex            = typename PowerBasisFE<Basis>::GlobalIndex;
    using FERequirementType      = FErequirements<Eigen::VectorXd>;
    using ResultRequirementsType = ResultRequirements<Eigen::VectorXd>;
    using LocalView              = typename Basis::LocalView;
    using GridView               = typename Basis::GridView;

    using Traits = TraitsFromLocalView<LocalView>;

    static constexpr int mydim = Traits::mydim;

    template <typename VolumeLoad = std::nullptr_t, typename NeumannBoundaryLoad = std::nullptr_t>
    requires(Std::is_pointer<VolumeLoad>and Std::is_pointer<NeumannBoundaryLoad>)
        LinearElastic(Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu,
                      VolumeLoad p_volumeLoad = nullptr, const BoundaryPatch<GridView>* neumannBoundary = nullptr,
                      NeumannBoundaryLoad p_neumannBoundaryLoad = nullptr)
        : BaseDisp(globalBasis, element),
          localView_{globalBasis.localView()},
          volumeLoad{Std::returnReferenceOrNulloptIfObjectIsNullPtr(p_volumeLoad)},
          neumannBoundaryLoad{Std::returnReferenceOrNulloptIfObjectIsNullPtr(p_neumannBoundaryLoad)},
          neumannBoundary_{neumannBoundary},
          emod_{emod},
          nu_{nu} {
      localView_.bind(element);
      auto& first_child = localView_.tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes     = fe.size();
      dispAtNodes.resize(fe.size());
      const int order = 2 * (localView_.tree().child(0).finiteElement().localBasis().order());
      localBasis      = Dune::CachedLocalBasis(localView_.tree().child(0).finiteElement().localBasis());
      localBasis.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order),
                      Dune::bindDerivatives(0, 1));

      assert(((not neumannBoundary_ and not neumannBoundaryLoad) or (neumannBoundary_ and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

  public:
    const auto& localView() const { return localView_; }

    auto getDisplacementFunction(const FERequirementType& par) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);

      for (auto i = 0U; i < dispAtNodes.size(); ++i)
        for (auto k2 = 0U; k2 < mydim; ++k2)
          dispAtNodes[i][k2] = d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];

      auto geo = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          localView_.element().geometry());
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

      const auto geo = localView_.element().geometry();
      double energy  = 0.0;
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const auto EVoigt = eps.evaluate(gpIndex, on(gridElement));

        energy += (0.5 * EVoigt.dot(C * EVoigt)) * geo.integrationElement(gp.position()) * gp.weight();
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
          const auto uVal                              = u.evaluate(gpIndex);
          Eigen::Vector<double, Traits::worlddim> fext = (*volumeLoad)(toEigen(gp.position()), lambda);
          energy -= uVal.dot(fext) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e. neumann boundary
      if (not neumannBoundary_ and not neumannBoundaryLoad) return energy;

      auto element = localView_.element();
      for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
        if (not neumannBoundary_->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, mydim - 1>::rule(intersection.type(), u.order());

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, mydim>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto uVal = u.evaluate(quadPos);

          // Value of the Neumann data at the current position
          auto neumannValue
              = (*neumannBoundaryLoad)(toEigen(intersection.geometry().global(curQuad.position())), lambda);

          energy -= neumannValue.dot(uVal) * curQuad.weight() * integrationElement;
        }
      }
      return energy;
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& K) const {
      const auto eps = getStrainFunction(par);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto C   = getMaterialTangent();
      const auto geo = localView_.element().geometry();

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

      auto cauchyStress = C * (epsVoigt).eval();

      typename ResultTypeMap<double>::ResultArray resultVector;
      if (req.isResultRequested(ResultType::cauchyStress)) {
        resultVector.resize(3, 1);
        resultVector = cauchyStress;
        result.insertOrAssignResult(ResultType::cauchyStress, resultVector);
      }
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& force) const {
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto eps     = getStrainFunction(par);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto C   = getMaterialTangent();
      const auto geo = localView_.element().geometry();

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
          Eigen::Vector<double, Traits::worlddim> fext = (*volumeLoad)(toEigen(gp.position()), lambda);
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<mydim>(mydim * i)
                -= udCi * fext * geo.integrationElement(gp.position()) * gp.weight();
          }
        }
      }

      // External forces, boundary forces, i.e. at the Neumann boundary
      if (not neumannBoundary_ and not neumannBoundaryLoad) return;

      const auto u = getDisplacementFunction(par);
      auto element = localView_.element();
      for (auto&& intersection : intersections(neumannBoundary_->gridView(), element)) {
        if (not neumannBoundary_->contains(intersection)) continue;

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
                = (*neumannBoundaryLoad)(toEigen(intersection.geometry().global(curQuad.position())), lambda);
            force.template segment<mydim>(mydim * i) -= udCi * neumannValue * curQuad.weight() * integrationElement;
          }
        }
      }
    }

    LocalView localView_;
    Dune::CachedLocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::optional<std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                                        const double&)>>
        volumeLoad;
    std::optional<std::function<Eigen::Vector<double, Traits::worlddim>(const Eigen::Vector<double, Traits::worlddim>&,
                                                                        const double&)>>
        neumannBoundaryLoad;
    const BoundaryPatch<GridView>* neumannBoundary_;
    mutable Dune::BlockVector<Dune::RealTuple<double, Traits::dimension>> dispAtNodes;
    double emod_;
    double nu_;
    size_t numberOfNodes{0};
  };

}  // namespace Ikarus
