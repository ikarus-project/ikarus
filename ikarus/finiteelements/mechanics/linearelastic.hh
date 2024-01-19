// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file linearelastic.hh
 * @brief Definition of the LinearElastic class for finite element mechanics computations.
 * @ingroup  mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS

#  include <iosfwd>
#  include <optional>
#  include <type_traits>

#  include <dune/fufem/boundarypatch.hh>
#  include <dune/geometry/quadraturerules.hh>
#  include <dune/localfefunctions/expressions/linearStrainsExpr.hh>
#  include <dune/localfefunctions/impl/standardLocalFunction.hh>
#  include <dune/localfefunctions/manifolds/realTuple.hh>

#  include <ikarus/finiteelements/febases/powerbasisfe.hh>
#  include <ikarus/finiteelements/fehelper.hh>
#  include <ikarus/finiteelements/ferequirements.hh>
#  include <ikarus/finiteelements/mechanics/loads.hh>
#  include <ikarus/finiteelements/mechanics/materials.hh>
#  include <ikarus/finiteelements/physicshelper.hh>
#  include <ikarus/utils/defaultfunctions.hh>
#  include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

  /**
   * @brief LinearElastic class represents a linear elastic finite element.
   *
   * @ingroup mechanics
   *
   * @tparam Basis_ The basis type for the finite element.
   * @tparam FERequirements_ The requirements for the finite element.
   * @tparam useEigenRef A boolean flag indicating whether to use Eigen references.
   */
  template <typename Basis_, typename FERequirements_ = FERequirements<>, bool useEigenRef = false>
  class LinearElastic : public PowerBasisFE<typename Basis_::FlatBasis>,
                        public Volume<LinearElastic<Basis_, FERequirements_, useEigenRef>,
                                      TraitsFromFE<Basis_, FERequirements_, useEigenRef>>,
                        public Traction<LinearElastic<Basis_, FERequirements_, useEigenRef>,
                                        TraitsFromFE<Basis_, FERequirements_, useEigenRef>> {
  public:
    using Traits                 = TraitsFromFE<Basis_, FERequirements_, useEigenRef>;
    using Basis                  = typename Traits::Basis;
    using FlatBasis              = typename Traits::FlatBasis;
    using FERequirementType      = typename Traits::FERequirementType;
    using LocalView              = typename Traits::LocalView;
    using Geometry               = typename Traits::Geometry;
    using GridView               = typename Traits::GridView;
    using Element                = typename Traits::Element;
    using ResultRequirementsType = typename Traits::ResultRequirementsType;
    using BaseDisp               = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using VolumeType             = Volume<LinearElastic<Basis_, FERequirements_, useEigenRef>, Traits>;
    using TractionType           = Traction<LinearElastic<Basis_, FERequirements_, useEigenRef>, Traits>;
    static constexpr int myDim   = Traits::mydim;
    using LocalBasisType         = decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis());

    /**
     * @brief Constructor for the LinearElastic class.
     *
     * @tparam VolumeLoad The type for the volume load function.
     * @tparam NeumannBoundaryLoad The type for the Neumann boundary load function.
     * @param globalBasis The global basis for the finite element.
     * @param element The element for which the finite element is constructed.
     * @param emod Young's modulus.
     * @param nu Poisson's ratio.
     * @param p_volumeLoad Volume load function (default is LoadDefault).
     * @param p_neumannBoundary Neumann boundary patch (default is nullptr).
     * @param p_neumannBoundaryLoad Neumann boundary load function (default is LoadDefault).
     */
    template <typename VolumeLoad = utils::LoadDefault, typename NeumannBoundaryLoad = utils::LoadDefault>
    LinearElastic(const Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu,
                  VolumeLoad p_volumeLoad = {}, const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
                  NeumannBoundaryLoad p_neumannBoundaryLoad = {})
        : BaseDisp(globalBasis.flat(), element),
          VolumeType(p_volumeLoad),
          TractionType(p_neumannBoundary, p_neumannBoundaryLoad),
          emod_{emod},
          nu_{nu} {
      this->localView().bind(element);
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      geo_              = std::make_shared<const Geometry>(this->localView().element().geometry());
      numberOfNodes_    = fe.size();
      order_            = 2 * (fe.localBasis().order());
      localBasis        = Dune::CachedLocalBasis(fe.localBasis());
      if constexpr (requires { this->localView().element().impl().getQuadratureRule(order_); })
        if (this->localView().element().impl().isTrimmed())
          localBasis.bind(this->localView().element().impl().getQuadratureRule(order_), Dune::bindDerivatives(0, 1));
        else
          localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order_),
                          Dune::bindDerivatives(0, 1));
      else
        localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order_),
                        Dune::bindDerivatives(0, 1));
    }
    /**
     * @brief Gets the displacement function for the given FERequirementType and optional displacement vector.
     *
     * @tparam ScalarType The scalar type for the displacement vector.
     * @param par The FERequirementType object.
     * @param dx Optional displacement vector.
     * @return The displacement function.
     */
    template <typename ScalarType = double>
    auto displacementFunction(const FERequirementType& par,
                              const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);
      auto disp     = Ikarus::FEHelper::localSolutionBlockVector<Traits>(d, this->localView(), dx);
      Dune::StandardLocalFunction uFunction(localBasis, disp, geo_);
      return uFunction;
    }
    /**
     * @brief Gets the strain function for the given FERequirementType and optional displacement vector.
     *
     * @tparam ScalarType The scalar type for the strain vector.
     * @param par The FERequirementType object.
     * @param dx Optional displacement vector.
     * @return The strain function.
     */
    template <class ScalarType = double>
    auto strainFunction(const FERequirementType& par,
                        const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      return Dune::linearStrains(displacementFunction(par, dx));
    }

    /**
     * @brief Gets the material tangent matrix for the linear elastic material.
     *
     * @return The material tangent matrix.
     */
    auto materialTangent() const {
      if constexpr (myDim == 2)
        return planeStressLinearElasticMaterialTangent(emod_, nu_);
      else if constexpr (myDim == 3)
        return linearElasticMaterialTangent3D(emod_, nu_);
    }

    /**
     * @brief Gets the material tangent function for the given FERequirementType.
     *
     * @param par The FERequirementType object.
     * @return The material tangent function.
     */
    auto materialTangentFunction([[maybe_unused]] const FERequirementType& par) const {
      return [&]([[maybe_unused]] auto gp) { return materialTangent(); };
    }

    const Geometry& geometry() const { return *geo_; }
    [[nodiscard]] size_t numberOfNodes() const { return numberOfNodes_; }
    [[nodiscard]] int order() const { return order_; }

    /**
     * @brief Calculates the scalar energy for the given FERequirementType.
     *
     * @param par The FERequirementType object.
     * @return The scalar energy.
     */
    inline double calculateScalar(const FERequirementType& par) const { return calculateScalarImpl<double>(par); }

    /**
     * @brief Calculates the vector force for the given FERequirementType.
     *
     * @param par The FERequirementType object.
     * @param force Vector to store the calculated force.
     */
    inline void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(par, force);
    }

    /**
     * @brief Calculates the matrix stiffness for the given FERequirementType.
     *
     * @param par The FERequirementType object.
     * @param K Matrix to store the calculated stiffness.
     */
    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> K) const {
      const auto eps = strainFunction(par);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto C = materialTangent();
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
        for (size_t i = 0; i < numberOfNodes_; ++i) {
          const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          for (size_t j = 0; j < numberOfNodes_; ++j) {
            const auto bopJ = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), on(gridElement));
            K.template block<myDim, myDim>(i * myDim, j * myDim) += bopI.transpose() * C * bopJ * intElement;
          }
        }
      }

      // Update due to displacement-dependent external loads, e.g., follower loads
      VolumeType::calculateMatrix(par, K);
      TractionType::calculateMatrix(par, K);
    }

    /**
     * @brief Calculates results at a specific local position.
     *
     * @param req The ResultRequirementsType object specifying the requested results.
     * @param local Local position vector.
     * @param result Map to store the calculated results.
     */
    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::Indices;
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto eps = strainFunction(req.getFERequirements());
      const auto C   = materialTangent();
      auto epsVoigt  = eps.evaluate(local, on(gridElement));

      auto linearStress = (C * epsVoigt).eval();

      if (req.isResultRequested(ResultType::linearStress))
        result.insertOrAssignResult(ResultType::linearStress, linearStress);
      else
        DUNE_THROW(Dune::NotImplemented, "The requested result type is NOT implemented.");
    }

  private:
    std::shared_ptr<const Geometry> geo_;
    Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisType>> localBasis;
    double emod_;
    double nu_;
    size_t numberOfNodes_{0};
    int order_{};

  protected:
    template <typename ScalarType>
    auto calculateScalarImpl(const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                                           = std::nullopt) const -> ScalarType {
      const auto uFunction = displacementFunction(par, dx);
      const auto eps       = strainFunction(par, dx);
      const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto C = materialTangent();

      ScalarType energy = 0.0;
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const auto EVoigt = eps.evaluate(gpIndex, on(gridElement));
        energy += (0.5 * EVoigt.dot(C * EVoigt)) * geo_->integrationElement(gp.position()) * gp.weight();
      }

      // External forces volume forces over the domain
      energy += VolumeType::calculateScalarImpl(par, dx);

      // line or surface loads, i.e., neumann boundary
      energy += TractionType::calculateScalarImpl(par, dx);
      return energy;
    }

    template <typename ScalarType>
    void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto eps = strainFunction(par, dx);
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      const auto C = materialTangent();

      // Internal forces
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
        auto stresses           = (C * eps.evaluate(gpIndex, on(gridElement))).eval();
        for (size_t i = 0; i < numberOfNodes_; ++i) {
          const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          force.template segment<myDim>(myDim * i) += bopI.transpose() * stresses * intElement;
        }
      }

      // External forces volume forces over the domain
      VolumeType::calculateVectorImpl(par, force, dx);

      // External forces, boundary forces, i.e., at the Neumann boundary
      TractionType::calculateVectorImpl(par, force, dx);
    }
  };
}  // namespace Ikarus

#else
#  error LinearElastic depends on dune-localfefunctions, which is not included
#endif
