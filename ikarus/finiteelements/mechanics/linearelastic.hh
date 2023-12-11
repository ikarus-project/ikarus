// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#if HAVE_DUNE_LOCALFEFUNCTIONS
#  include <iosfwd>
#  include <optional>
#  include <type_traits>

#  include <dune/common/classname.hh>
#  include <dune/fufem/boundarypatch.hh>
#  include <dune/geometry/quadraturerules.hh>
#  include <dune/geometry/type.hh>
#  include <dune/localfefunctions/expressions/linearStrainsExpr.hh>
#  include <dune/localfefunctions/impl/standardLocalFunction.hh>
#  include <dune/localfefunctions/manifolds/realTuple.hh>

#  include <autodiff/forward/dual.hpp>
#  include <autodiff/forward/dual/eigen.hpp>

#  include <ikarus/finiteelements/febases/powerbasisfe.hh>
#  include <ikarus/finiteelements/ferequirements.hh>
#  include <ikarus/finiteelements/fetraits.hh>
#  include <ikarus/finiteelements/mechanics/materials.hh>
#  include <ikarus/finiteelements/physicshelper.hh>
#  include <ikarus/utils/defaultfunctions.hh>
#  include <ikarus/utils/eigendunetransformations.hh>
#  include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

  template <typename Basis_, typename FERequirements_ = FErequirements<>, bool useEigenRef = false>
  class LinearElastic : public PowerBasisFE<typename Basis_::FlatBasis> {
  public:
    using Basis                  = Basis_;
    using FlatBasis              = typename Basis::FlatBasis;
    using BaseDisp               = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType      = FERequirements_;
    using ResultRequirementsType = ResultRequirements<FERequirementType>;
    using LocalView              = typename FlatBasis::LocalView;
    using Element                = typename LocalView::Element;
    using GridView               = typename FlatBasis::GridView;

    using Traits = TraitsFromLocalView<LocalView, useEigenRef>;

    static constexpr int myDim = Traits::mydim;

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
      order      = 2 * (this->localView().tree().child(0).finiteElement().localBasis().order());
      localBasis = Dune::CachedLocalBasis(this->localView().tree().child(0).finiteElement().localBasis());
      if constexpr (requires { this->localView().element().impl().getQuadratureRule(order); })
        if (this->localView().element().impl().isTrimmed())
          localBasis.bind(this->localView().element().impl().getQuadratureRule(order), Dune::bindDerivatives(0, 1));
        else
          localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                          Dune::bindDerivatives(0, 1));
      else
        localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                        Dune::bindDerivatives(0, 1));

      if constexpr (!std::is_same_v<VolumeLoad, LoadDefault>) volumeLoad = p_volumeLoad;
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, LoadDefault>) neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

  public:
    template <typename ScalarType>
    auto getDisplacementFunction(const FERequirementType& par,
                                 const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);

      Dune::BlockVector<Dune::RealTuple<ScalarType, Traits::dimension>> disp(dispAtNodes.size());
      if (dx) {
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < myDim; ++k2)
            disp[i][k2] = dx.value()[i * myDim + k2]
                          + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      } else
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < myDim; ++k2)
            disp[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];

      auto geo = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          this->localView().element().geometry());
      Dune::StandardLocalFunction uFunction(localBasis, disp, geo);

      return uFunction;
    }

    template <class ScalarType = double>
    auto getStrainFunction(const FERequirementType& par,
                           const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      return linearStrains(getDisplacementFunction(par, dx));
    }

    auto getMaterialTangent() const {
      if constexpr (myDim == 2)
        return planeStressLinearElasticMaterialTangent(emod_, nu_);
      else if constexpr (myDim == 3)
        return linearElasticMaterialTangent3D(emod_, nu_);
    }

    auto getMaterialTangentFunction([[maybe_unused]] const FERequirementType& par) const {
      return [&]([[maybe_unused]] auto gp) { return getMaterialTangent(); };
    }

    inline double calculateScalar(const FERequirementType& par) const { return calculateScalarImpl<double>(par); }

    inline void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(par, force);
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> K) const {
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
            K.template block<myDim, myDim>(i * myDim, j * myDim) += bopI.transpose() * C * bopJ * intElement;
          }
        }
      }
    }

    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::Indices;
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto eps = getStrainFunction(req.getFERequirements());
      const auto C   = getMaterialTangent();
      auto epsVoigt  = eps.evaluate(local, on(gridElement));

      auto linearStress = (C * epsVoigt).eval();

      if (req.isResultRequested(ResultType::linearStress))
        result.insertOrAssignResult(ResultType::linearStress, linearStress);
      else
        DUNE_THROW(Dune::NotImplemented, "The requested result type is NOT implemented.");
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
    int order{};

  protected:
    template <typename ScalarType>
    auto calculateScalarImpl(const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                                           = std::nullopt) const -> ScalarType {
      const auto u       = getDisplacementFunction(par, dx);
      const auto eps     = getStrainFunction(par, dx);
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto C = getMaterialTangent();

      const auto geo    = this->localView().element().geometry();
      ScalarType energy = 0.0;
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const auto EVoigt = eps.evaluate(gpIndex, on(gridElement));

        energy += (0.5 * EVoigt.dot(C * EVoigt)) * geo.integrationElement(gp.position()) * gp.weight();
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
          const auto uVal                              = u.evaluate(gpIndex);
          Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigen(geo.global(gp.position())), lambda);
          energy -= uVal.dot(fext) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e., neumann boundary
      if (not neumannBoundary and not neumannBoundaryLoad) return energy;

      auto element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, myDim - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());

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

    template <typename ScalarType>
    void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto eps     = getStrainFunction(par, dx);
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
          force.template segment<myDim>(myDim * i) += bopI.transpose() * stresses * intElement;
        }
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        const auto u = getDisplacementFunction(par, dx);
        for (const auto& [gpIndex, gp] : u.viewOverIntegrationPoints()) {
          Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(toEigen(geo.global(gp.position())), lambda);
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<myDim>(myDim * i)
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
          const Dune::FieldVector<double, myDim>& quadPos = intersection.geometryInInside().global(curQuad.position());

          const double integrationElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function wrt the i-th coef
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = u.evaluateDerivative(quadPos, wrt(coeff(i)));

            // Value of the Neumann data at the current position
            auto neumannValue
                = neumannBoundaryLoad(toEigen(intersection.geometry().global(curQuad.position())), lambda);
            force.template segment<myDim>(myDim * i) -= udCi * neumannValue * curQuad.weight() * integrationElement;
          }
        }
      }
    }
  };
}  // namespace Ikarus

#else
#  error LinearElastic depends on dune-localfefunctions, which is not included
#endif
