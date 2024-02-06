// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearelastic.hh
 * \brief Definition of the NonLinearElastic class for finite element mechanics computations.
 * \ingroup  mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/fufem/boundarypatch.hh>
  #include <dune/geometry/quadraturerules.hh>
  #include <dune/geometry/type.hh>
  #include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
  #include <dune/localfefunctions/expressions/greenLagrangeStrains.hh>
  #include <dune/localfefunctions/impl/standardLocalFunction.hh>
  #include <dune/localfefunctions/manifolds/realTuple.hh>

  #include <ikarus/finiteelements/febases/powerbasisfe.hh>
  #include <ikarus/finiteelements/fehelper.hh>
  #include <ikarus/finiteelements/ferequirements.hh>
  #include <ikarus/finiteelements/fetraits.hh>
  #include <ikarus/finiteelements/mechanics/loads.hh>
  #include <ikarus/finiteelements/mechanics/materials/tags.hh>
  #include <ikarus/finiteelements/physicshelper.hh>
  #include <ikarus/utils/defaultfunctions.hh>
  #include <ikarus/utils/eigendunetransformations.hh>
  #include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

/**
 * \brief NonLinearElastic class represents a non-linear elastic finite element.
 *
 * \ingroup mechanics
 *
 * \tparam B The basis type for the finite element.
 * \tparam MAT The material type for the finite element.
 * \tparam FER The requirements for the finite element.
 * \tparam useEigenRef A boolean flag indicating whether to use Eigen references.
 */
template <typename B, typename MAT, typename FER = FERequirements<>, bool useEigenRef = false>
class NonLinearElastic : public PowerBasisFE<B>,
                         public Volume<NonLinearElastic<B, MAT, FER, useEigenRef>, FETraits<B, FER, useEigenRef>>,
                         public Traction<NonLinearElastic<B, MAT, FER, useEigenRef>, FETraits<B, FER, useEigenRef>>
{
public:
  using Traits            = FETraits<B, FER, useEigenRef>;
  using Basis             = typename Traits::Basis;
  using FlatBasis         = typename Traits::FlatBasis;
  using FERequirementType = typename Traits::FERequirementType;
  using LocalView         = typename Traits::LocalView;
  using Geometry          = typename Traits::Geometry;
  using GridView          = typename Traits::GridView;
  using Element           = typename Traits::Element;
  using BasePowerFE       = PowerBasisFE<Basis>; // Handles globalIndices function
  using Material          = MAT;
  using VolumeType        = Volume<NonLinearElastic, Traits>;
  using TractionType      = Traction<NonLinearElastic, Traits>;
  using LocalBasisType    = decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis());

  static constexpr int myDim       = Traits::mydim;
  static constexpr auto strainType = StrainTags::greenLagrangian;

  /**
   * \brief Constructor for the NonLinearElastic class.
   *
   * \tparam VolumeLoad The type for the volume load function.
   * \tparam NeumannBoundaryLoad The type for the Neumann boundary load function.
   * \param globalBasis The global basis for the finite element.
   * \param element The element for which the finite element is constructed.
   * \param mat The material for the non-linear elastic element.
   * \param volumeLoad Volume load function (default is LoadDefault).
   * \param neumannBoundary Neumann boundary patch (default is nullptr).
   * \param neumannBoundaryLoad Neumann boundary load function (default is LoadDefault).
   */
  template <typename VolumeLoad = utils::LoadDefault, typename NeumannBoundaryLoad = utils::LoadDefault>
  NonLinearElastic(const Basis& globalBasis, const typename LocalView::Element& element, const Material& mat,
                   VolumeLoad volumeLoad = {}, const BoundaryPatch<GridView>* neumannBoundary = nullptr,
                   NeumannBoundaryLoad neumannBoundaryLoad = {})
      : BasePowerFE(globalBasis, element),
        VolumeType(volumeLoad),
        TractionType(neumannBoundary, neumannBoundaryLoad),
        mat_{mat} {
    this->localView().bind(element);
    auto& firstChild = this->localView().tree().child(0);
    const auto& fe   = firstChild.finiteElement();
    geo_             = std::make_shared<const Geometry>(this->localView().element().geometry());
    numberOfNodes_   = fe.size();
    order_           = 2 * (fe.localBasis().order());
    localBasis_      = Dune::CachedLocalBasis(fe.localBasis());
    if constexpr (requires { this->localView().element().impl().getQuadratureRule(order_); })
      if (this->localView().element().impl().isTrimmed())
        localBasis_.bind(this->localView().element().impl().getQuadratureRule(order_), Dune::bindDerivatives(0, 1));
      else
        localBasis_.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order_),
                         Dune::bindDerivatives(0, 1));
    else
      localBasis_.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order_),
                       Dune::bindDerivatives(0, 1));
  }

  /**
   * \brief Get the displacement function for the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the displacement function.
   * \param par The FERequirementType object.
   * \param dx Optional displacement vector.
   * \return A StandardLocalFunction representing the displacement function.
   */
  template <typename ScalarType = double>
  auto displacementFunction(const FERequirementType& par,
                            const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
    const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);
    auto disp     = Ikarus::FEHelper::localSolutionBlockVector<Traits>(d, this->localView(), dx);
    Dune::StandardLocalFunction uFunction(localBasis_, disp, geo_);
    return uFunction;
  }

  /**
   * \brief The strain function for the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the strain function.
   * \param par The FERequirementType object.
   * \param dx Optional displacement vector.
   * \return The strain function calculated using greenLagrangeStrains.
   */
  template <typename ScalarType = double>
  inline auto strainFunction(const FERequirementType& par,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
    return Dune::greenLagrangeStrains(displacementFunction(par, dx));
  }

  /**
   * \brief Get the material tangent for the given strain.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \tparam voigt Flag indicating whether to use Voigt notation.
   * \param strain The strain vector.
   * \return The material tangent calculated using the material's tangentModuli function.
   */
  template <typename ScalarType, int strainDim, bool voigt = true>
  auto materialTangent(const Eigen::Vector<ScalarType, strainDim>& strain) const {
    if constexpr (std::is_same_v<ScalarType, double>)
      return mat_.template tangentModuli<strainType, voigt>(strain);
    else {
      decltype(auto) matAD = mat_.template rebind<ScalarType>();
      return matAD.template tangentModuli<strainType, voigt>(strain);
    }
  }

  /**
   * \brief Get the internal energy for the given strain.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \param strain The strain vector.
   * \return The internal energy calculated using the material's storedEnergy function.
   */
  template <typename ScalarType, int strainDim>
  auto getInternalEnergy(const Eigen::Vector<ScalarType, strainDim>& strain) const {
    if constexpr (std::is_same_v<ScalarType, double>)
      return mat_.template storedEnergy<strainType>(strain);
    else {
      decltype(auto) matAD = mat_.template rebind<ScalarType>();
      return matAD.template storedEnergy<strainType>(strain);
    }
  }

  /**
   * \brief Get the stress for the given strain.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \tparam voigt A boolean indicating whether to use the Voigt notation for stress.
   * \param strain The strain vector.
   * \return The stress vector calculated using the material's stresses function.
   */
  template <typename ScalarType, int strainDim, bool voigt = true>
  auto getStress(const Eigen::Vector<ScalarType, strainDim>& strain) const {
    if constexpr (std::is_same_v<ScalarType, double>)
      return mat_.template stresses<strainType, voigt>(strain);
    else {
      decltype(auto) matAD = mat_.template rebind<ScalarType>();
      return matAD.template stresses<strainType, voigt>(strain);
    }
  }

  const Geometry& geometry() const { return *geo_; }
  [[nodiscard]] size_t numberOfNodes() const { return numberOfNodes_; }
  [[nodiscard]] int order() const { return order_; }

  /**
   * \brief Calculate the scalar value associated with the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param par The FERequirementType object specifying the requirements for the calculation.
   * \return The calculated scalar value.
   */
  double calculateScalar(const FERequirementType& par) const { return calculateScalarImpl<double>(par); }

  /**
   * \brief Calculate the vector associated with the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param par The FERequirementType object specifying the requirements for the calculation.
   * \param force The vector to store the calculated result.
   */
  void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> force) const {
    calculateVectorImpl<double>(par, force);
  }

  /**
   * \brief Calculate the matrix associated with the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param par The FERequirementType object specifying the requirements for the calculation.
   * \param K The matrix to store the calculated result.
   */
  void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> K) const {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto eps = strainFunction(par);
    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
      const auto EVoigt       = (eps.evaluate(gpIndex, on(gridElement))).eval();
      const auto C            = materialTangent(EVoigt);
      const auto stresses     = getStress(EVoigt);
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        for (size_t j = 0; j < numberOfNodes_; ++j) {
          const auto bopJ = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), on(gridElement));
          const auto kgIJ = eps.evaluateDerivative(gpIndex, wrt(coeff(i, j)), along(stresses), on(gridElement));
          K.template block<myDim, myDim>(i * myDim, j * myDim) += (bopI.transpose() * C * bopJ + kgIJ) * intElement;
        }
      }
    }

    // Update due to displacement-dependent external loads, e.g., follower loads
    VolumeType::calculateMatrix(par, K);
    TractionType::calculateMatrix(par, K);
  }

  /**
   * \brief Calculates a requested result at a specific local position.
   *
   * \param req The FERequirementType object holding the global solution.
   * \param local Local position vector.
   * \return calculated result
   *
   * \tparam resType The type representing the requested result.
   */
  template <ResultType resType>
  auto calculateAt(const FERequirementType& req, const Dune::FieldVector<double, Traits::mydim>& local) const {
    static_assert(resType == ResultType::PK2Stress, "The requested result type is NOT implemented.");

    using namespace Dune::DerivativeDirections;
    if constexpr (resType == ResultType::PK2Stress) {
      const auto uFunction = displacementFunction(req);
      const auto H         = uFunction.evaluateDerivative(local, Dune::wrt(spatialAll), Dune::on(gridElement));
      const auto E         = (0.5 * (H.transpose() + H + H.transpose() * H)).eval();
      const auto EVoigt    = toVoigt(E);

      return mat_.template stresses<StrainTags::greenLagrangian>(EVoigt);
    }
  }

private:
  std::shared_ptr<const Geometry> geo_;
  Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisType>> localBasis_;
  Material mat_;
  size_t numberOfNodes_{0};
  int order_{};

protected:
  template <typename ScalarType>
  auto calculateScalarImpl(const FERequirementType& par,
                           const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const
      -> ScalarType {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto uFunction = displacementFunction(par, dx);
    const auto eps       = strainFunction(par, dx);
    const auto& lambda   = par.getParameter(Ikarus::FEParameter::loadfactor);
    ScalarType energy    = 0.0;

    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const auto EVoigt         = (eps.evaluate(gpIndex, on(gridElement))).eval();
      const auto internalEnergy = getInternalEnergy(EVoigt);
      energy += internalEnergy * geo_->integrationElement(gp.position()) * gp.weight();
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
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto eps = strainFunction(par, dx);

    // Internal forces
    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
      const auto EVoigt       = (eps.evaluate(gpIndex, on(gridElement))).eval();
      const auto stresses     = getStress(EVoigt);
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
} // namespace Ikarus

#else
  #error NonLinearElastic depends on dune-localfefunctions, which is not included
#endif
