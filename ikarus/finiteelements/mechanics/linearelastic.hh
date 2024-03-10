// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file linearelastic.hh
 * \brief Definition of the LinearElastic class for finite element mechanics computations.
 * \ingroup  mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS

  #include <optional>
  #include <type_traits>

  #include <dune/geometry/quadraturerules.hh>
  #include <dune/localfefunctions/expressions/linearStrainsExpr.hh>
  #include <dune/localfefunctions/impl/standardLocalFunction.hh>

  #include <ikarus/finiteelements/fehelper.hh>
  #include <ikarus/finiteelements/ferequirements.hh>
  #include <ikarus/finiteelements/mechanics/materials.hh>
  #include <ikarus/finiteelements/physicshelper.hh>

namespace Ikarus {

template <typename PreFE, typename FE>
class LinearElastic;

/**
 * \brief A PreFE struct for linear elastic elements.
 */
struct LinearElasticPre
{
  YoungsModulusAndPoissonsRatio material;

  template <typename PreFE, typename FE>
  using Skill = LinearElastic<PreFE, FE>;
};

/**
 * \brief LinearElastic class represents a linear elastic finite element.
 *
 * \ingroup mechanics
 *
 * \tparam PreFE The type of the  pre finite element.
 * \tparam FE The type of the finite element.
 */
template <typename PreFE, typename FE>
class LinearElastic
{
public:
  using Traits            = PreFE::Traits;
  using BasisHandler      = typename Traits::BasisHandler;
  using FlatBasis         = typename Traits::FlatBasis;
  using FERequirementType = typename Traits::FERequirementType;
  using LocalView         = typename Traits::LocalView;
  using Geometry          = typename Traits::Geometry;
  using GridView          = typename Traits::GridView;
  using Element           = typename Traits::Element;
  using Pre               = LinearElasticPre;

  static constexpr int myDim = Traits::mydim;
  using LocalBasisType       = decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis());

  /**
   * \brief Constructor for the LinearElastic class.
   * \param pre The pre fe
   */
  explicit LinearElastic(const Pre& pre)
      : mat_{pre.material} {}

protected:
  /**
   * \brief A helper function to bind the local view to the element.
   */
  void bindImpl() {
    const auto& localView = underlying().localView();
    const auto& element   = localView.element();
    auto& firstChild      = localView.tree().child(0);
    const auto& fe        = firstChild.finiteElement();
    geo_                  = std::make_shared<const Geometry>(element.geometry());
    numberOfNodes_        = fe.size();
    order_                = 2 * (fe.localBasis().order());
    localBasis_           = Dune::CachedLocalBasis(fe.localBasis());
    if constexpr (requires { element.impl().getQuadratureRule(order_); })
      if (element.impl().isTrimmed())
        localBasis_.bind(element.impl().getQuadratureRule(order_), Dune::bindDerivatives(0, 1));
      else
        localBasis_.bind(Dune::QuadratureRules<double, myDim>::rule(element.type(), order_),
                         Dune::bindDerivatives(0, 1));
    else
      localBasis_.bind(Dune::QuadratureRules<double, myDim>::rule(element.type(), order_), Dune::bindDerivatives(0, 1));
  }

public:
  /**
   * \brief Gets the displacement function for the given FERequirementType and optional displacement vector.
   *
   * \tparam ScalarType The scalar type for the displacement vector.
   * \param par The FERequirementType object.
   * \param dx Optional displacement vector.
   * \return The displacement function.
   */
  template <typename ScalarType = double>
  auto displacementFunction(
      const FERequirementType& par,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    const auto& d = par.getGlobalSolution(Ikarus::FESolutions::displacement);
    auto disp     = Ikarus::FEHelper::localSolutionBlockVector<Traits>(d, underlying().localView(), dx);
    Dune::StandardLocalFunction uFunction(localBasis_, disp, geo_);
    return uFunction;
  }
  /**
   * \brief Gets the strain function for the given FERequirementType and optional displacement vector.
   *
   * \tparam ScalarType The scalar type for the strain vector.
   * \param par The FERequirementType object.
   * \param dx Optional displacement vector.
   * \return The strain function.
   */
  template <class ScalarType = double>
  auto strainFunction(
      const FERequirementType& par,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    return Dune::linearStrains(displacementFunction(par, dx));
  }

  /**
   * \brief Gets the material tangent matrix for the linear elastic material.
   *
   * \return The material tangent matrix.
   */
  auto materialTangent() const {
    if constexpr (myDim == 2)
      return planeStressLinearElasticMaterialTangent(mat_.emodul, mat_.nu);
    else if constexpr (myDim == 3)
      return linearElasticMaterialTangent3D(mat_.emodul, mat_.nu);
  }

  /**
   * \brief Gets the material tangent function for the given FERequirementType.
   *
   * \param par The FERequirementType object.
   * \return The material tangent function.
   */
  auto materialTangentFunction([[maybe_unused]] const FERequirementType& par) const {
    return [&]([[maybe_unused]] auto gp) { return materialTangent(); };
  }

  const Geometry& geometry() const { return *geo_; }
  [[nodiscard]] size_t numberOfNodes() const { return numberOfNodes_; }
  [[nodiscard]] int order() const { return order_; }

  using SupportedResultTypes = std::tuple<decltype(makeRT<ResultTypes::linearStress>())>;

  template <template <typename, int, int> class RT>
  static consteval bool canProvideResultType() {
    return isSameResultType<RT, ResultTypes::linearStress>;
  }

public:
  /**
   * \brief Calculates a requested result at a specific local position.
   *
   * \param req The FERequirementType object holding the global solution.
   * \param local Local position vector.
   * \tparam RT The requested result type
   * \return calculated result
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  requires(canProvideResultType<RT>())
  auto calculateAtImpl(const FERequirementType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<1>) const {
    using RTWrapper = ResultWrapper<RT<typename Traits::ctype, myDim, Traits::worlddim>, ResultShape::Vector>;
    if constexpr (isSameResultType<RT, ResultTypes::linearStress>) {
      const auto eps = strainFunction(req);
      const auto C   = materialTangent();
      auto epsVoigt  = eps.evaluate(local, Dune::on(Dune::DerivativeDirections::gridElement));

      return RTWrapper{(C * epsVoigt).eval()};
    }
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }

  std::shared_ptr<const Geometry> geo_;
  Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisType>> localBasis_;
  YoungsModulusAndPoissonsRatio mat_;
  size_t numberOfNodes_{0};
  int order_{};

protected:
  template <typename ScalarType>
  void calculateMatrixImpl(
      const FERequirementType& par, typename Traits::template MatrixType<> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    const auto eps = strainFunction(par, dx);
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
  }

  template <typename ScalarType>
  auto calculateScalarImpl(const FERequirementType& par,
                           const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx =
                               std::nullopt) const -> ScalarType {
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
    return energy;
  }

  template <typename ScalarType>
  void calculateVectorImpl(
      const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
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
  }
};

/**
 * \brief A helper function to create a linear elastic pre finite element.
 * \param mat Material parameters for the linear elastic element.
 * \return A linear elastic pre finite element.
 */
auto linearElastic(const YoungsModulusAndPoissonsRatio& mat) {
  LinearElasticPre pre(mat);

  return pre;
}
} // namespace Ikarus

#else
  #error LinearElastic depends on dune-localfefunctions, which is not included
#endif
