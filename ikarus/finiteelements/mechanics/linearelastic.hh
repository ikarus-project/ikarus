// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
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
  #include <ikarus/finiteelements/feresulttypes.hh>
  #include <ikarus/finiteelements/mechanics/materials.hh>
  #include <ikarus/finiteelements/mechanics/materials/tags.hh>
  #include <ikarus/finiteelements/mechanics/strainenhancements/easfunctions.hh>
  #include <ikarus/finiteelements/physicshelper.hh>

namespace Ikarus {

template <typename PreFE, typename FE, typename PRE>
class LinearElastic;

/**
 * \brief A PreFE struct for linear elastic elements.
 * \tparam MAT Type of the material.
 */
template <Concepts::GeometricallyLinearMaterial MAT>
struct LinearElasticPre
{
  using Material = MAT;
  MAT material;

  template <typename PreFE, typename FE>
  using Skill = LinearElastic<PreFE, FE, LinearElasticPre>;
};

/**
 * \brief LinearElastic class represents a linear elastic finite element.
 *
 * \ingroup mechanics
 *
 * \tparam PreFE The type of the  pre finite element.
 * \tparam FE The type of the finite element.
 */
template <typename PreFE, typename FE, typename PRE>
class LinearElastic : public ResultTypeBase<ResultTypes::linearStress, ResultTypes::linearStressFull>
{
public:
  using Traits       = PreFE::Traits;
  using BasisHandler = typename Traits::BasisHandler;
  using FlatBasis    = typename Traits::FlatBasis;
  using Requirement  = FERequirements<FESolutions::displacement, FEParameter::loadfactor>;
  using LocalView    = typename Traits::LocalView;
  using Geometry     = typename Traits::Geometry;
  using GridView     = typename Traits::GridView;
  using Element      = typename Traits::Element;
  using Material     = PRE::Material;
  using Pre          = PRE;

  using LocalBasisType = decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis());

  template <typename ST>
  using VectorXOptRef = std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>;

  static constexpr int myDim       = Traits::mydim;
  static constexpr int strainDim   = myDim * (myDim + 1) / 2;
  static constexpr auto strainType = StrainTags::linear;
  static constexpr auto stressType = StressTags::linear;

  template <template <typename, int, int> class RT>
  using RTWrapperType = ResultWrapper<RT<typename Traits::ctype, myDim, Traits::worlddim>, ResultShape::Vector>;

  template <typename ST = double>
  using StrainType = Eigen::Vector<ST, strainDim>;

  template <typename ST = double>
  using BopType = Eigen::Matrix<ST, strainDim, myDim>;

  template <typename ST = double>
  using KgType = Eigen::Matrix<ST, myDim, myDim>;

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
   * \brief Gets the displacement function for the given Requirement and optional displacement vector.
   *
   * \tparam ScalarType The scalar type for the displacement vector.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return The displacement function.
   */
  template <typename ScalarType = double>
  auto displacementFunction(const Requirement& par, const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    const auto& d = par.globalSolution();
    auto disp     = Ikarus::FEHelper::localSolutionBlockVector<Traits>(d, underlying().localView(), dx);
    Dune::StandardLocalFunction uFunction(localBasis_, disp, geo_);
    return uFunction;
  }
  /**
   * \brief Gets the strain function for the given Requirement and optional di
   splacement vector.
   *
   * \tparam ScalarType The scalar type for the strain vector.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return The strain function.
   */
  template <class ScalarType = double>
  auto strainFunction(const Requirement& par, const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    return Dune::linearStrains(displacementFunction(par, dx));
  }

  /**
   * \brief Get the internal energy for the given strain.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \param strain The strain vector in Voigt notation.
   * \return The internal energy calculated using the material's storedEnergy function.
   */
  template <typename ScalarType, int strainDim>
  auto internalEnergy(const Eigen::Vector<ScalarType, strainDim>& strain) const {
    const auto C = materialTangent<ScalarType, strainDim>(strain);
    return (0.5 * strain.dot(C * strain));
  }

  /**
   * \brief Get the stress for the given strain.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \tparam voigt A boolean indicating whether to use the Voigt notation for stress.
   * \param strain The strain vector in Voigt notation.
   * \return The stress vector calculated using the material's stresses function.
   */
  template <typename ScalarType, int strainDim, bool voigt = true>
  auto stress(const Eigen::Vector<ScalarType, strainDim>& strain) const {
    const auto C = materialTangent<ScalarType, strainDim, voigt>(strain);
    return (C * strain).eval();
  }

  /**
   * \brief Get the linear elastic material tangent for the given strain for the given Requirement.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \tparam strainDim The dimension of the strain vector.
   * \tparam voigt Flag indicating whether to use Voigt notation.
   * \param strain The strain vector.
   * \return The material tangent calculated using the material's tangentModuli function.
   */
  template <typename ScalarType, int strainDim, bool voigt = true>
  auto materialTangent(
      const Eigen::Vector<ScalarType, strainDim>& strain = Eigen::Vector<ScalarType, strainDim>::Zero()) const {
    // Since that material is independent of the strains, a zero strain is passed here
    return material<ScalarType>().template tangentModuli<strainType, voigt>(strain);
  }

  const Geometry& geometry() const { return *geo_; }
  [[nodiscard]] size_t numberOfNodes() const { return numberOfNodes_; }
  [[nodiscard]] int order() const { return order_; }
  const Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisType>>& localBasis() const { return localBasis_; }

  template <typename ScalarType = double>
  decltype(auto) material() const {
    if constexpr (Concepts::AutodiffScalar<ScalarType>)
      return mat_.template rebind<ScalarType>();
    else
      return mat_;
  }

public:
  /**
   * \brief Get a lambda function that evaluates the requested result type for a given strain (in Voigt notation).
   * \tparam RT The type representing the requested result.
   * \return A lambda function that evaluates the requested result type for a given strain (in Voigt notation).
   */
  template <template <typename, int, int> class RT>
  requires(LinearElastic::template supportsResultType<RT>())
  auto resultFunction() const {
    return [&](const Eigen::Vector<double, strainDim>& strainInVoigt) {
      if constexpr (isSameResultType<RT, ResultTypes::linearStress> or
                    isSameResultType<RT, ResultTypes::linearStressFull>) {
        decltype(auto) mat = [&]() {
          if constexpr (isSameResultType<RT, ResultTypes::linearStressFull> and requires { mat_.underlying(); })
            return mat_.underlying();
          else
            return mat_;
        }();
        return RTWrapperType<RT>{mat.template stresses<StrainTags::linear>(enlargeIfReduced<Material>(strainInVoigt))};
      }
    };
  }

  /**
   * \brief Calculates a requested result at a specific local position.
   * \param req The Requirement object holding the global solution.
   * \param local Local position vector.
   * \return Calculated result.
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  requires(LinearElastic::template supportsResultType<RT>())
  auto calculateAtImpl(const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<1>) const {
    if constexpr (FE::isMixed())
      return RTWrapperType<RT>{};
    if constexpr (isSameResultType<RT, ResultTypes::linearStress> or
                  isSameResultType<RT, ResultTypes::linearStressFull>) {
      const auto rFunction = resultFunction<RT>();
      const auto eps       = strainFunction(req);
      auto epsVoigt        = eps.evaluate(local, Dune::on(Dune::DerivativeDirections::gridElement));
      return rFunction(epsVoigt);
    }
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }

  std::shared_ptr<const Geometry> geo_;
  Dune::CachedLocalBasis<std::remove_cvref_t<LocalBasisType>> localBasis_;
  Material mat_;
  size_t numberOfNodes_{0};
  int order_{};

public:
  /**
   * \brief Get a lambda function that evaluates the geometric part of the stiffness matrix (Kg) for a given integration
   * point and its index.
   * \details This is here due to symmetry to the `NonLinearElastic` class
   *
   * \tparam ST The scalar type for the material and strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \param K The matrix to store the calculated result.
   * \return A lambda function that evaluates the stiffness matrix for a given strain, integration point and its index.
   */
  template <typename ST>
  auto geometricStiffnessMatrixFunction(const Requirement& par, typename Traits::template MatrixType<ST>& K,
                                        const VectorXOptRef<ST>& dx = std::nullopt) const {
    return [&](const FE::template KgType<ST>& kgIJ, const int I, const int J, const auto& gp) {
      const auto geo          = underlying().localView().element().geometry();
      const double intElement = geo.integrationElement(gp.position()) * gp.weight();
      K.template block<FE::myDim, FE::myDim>(I * FE::myDim, J * FE::myDim) += kgIJ * intElement;
    };
  }
  /**
   * \brief Get a lambda function that evaluates the elastic stiffness matrix for a given strain, integration point and
   * its index.
   *
   * \tparam ST The scalar type for the material and strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \param K The matrix to store the calculated result.
   * \return A lambda function that evaluates the stiffness matrix for a given strain, integration point and its index.
   */
  template <typename ST>
  auto materialStiffnessMatrixFunction(const Requirement& par, typename Traits::template MatrixType<ST>& K,
                                       const VectorXOptRef<ST>& dx = std::nullopt) const {
    return [&](const StrainType<ST>& strain, const BopType<ST>& bopI, const BopType<ST>& bopJ, const int I, const int J,
               const auto& gp) {
      const auto C            = materialTangent(strain);
      const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
      K.template block<myDim, myDim>(I * myDim, J * myDim) += (bopI.transpose() * C * bopJ) * intElement;
    };
  }

  /**
   * \brief Get a lambda function that evaluates the internal force vector for a given strain, integration point and its
   * index.
   *
   * \tparam ST The scalar type for the material and strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \param force The vector to store the calculated result.
   * \return A lambda function that evaluates the intenral force vector for a given strain, integration point and its
   * index.
   */
  template <typename ST>
  auto internalForcesFunction(const Requirement& par, typename Traits::template VectorType<ST>& force,
                              const VectorXOptRef<ST>& dx = std::nullopt) const {
    return [&](const StrainType<ST>& stresses, const BopType<ST>& bopI, const int I, const auto& gp) {
      const double intElement = geo_->integrationElement(gp.position()) * gp.weight();
      force.template segment<myDim>(myDim * I) += bopI.transpose() * stresses * intElement;
    };
  }

  /**
   * \brief Get a lambda function that evaluates the internal energy at a given integration point and its index.
   *
   * \tparam ScalarType The scalar type for the material and strain.
   * \param par The Requirement object.
   * \param dx Optional displacement vector.
   * \return A lambda function that returns the intenral energy at a given integration point and its index.
   */
  template <typename ScalarType>
  auto energyFunction(const Requirement& par, const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    return [&]() -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      ScalarType energy = 0.0;
      const auto eps    = strainFunction(par, dx);
      for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
        const auto epsVoigt = eps.evaluate(gpIndex, on(gridElement));
        energy += internalEnergy(epsVoigt) * geo_->integrationElement(gp.position()) * gp.weight();
      }
      return energy;
    };
  }

protected:
  template <typename ScalarType>
  void calculateMatrixImpl(const Requirement& par, const MatrixAffordance& affordance,
                           typename Traits::template MatrixType<> K,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    if constexpr (FE::isMixed()) {
      return;
    }
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto eps       = strainFunction(par, dx);
    const auto kFunction = materialStiffnessMatrixFunction<ScalarType>(par, K, dx);

    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const auto epsVoigt = eps.evaluate(gpIndex, on(gridElement));
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        for (size_t j = 0; j < numberOfNodes_; ++j) {
          const auto bopJ = eps.evaluateDerivative(gpIndex, wrt(coeff(j)), on(gridElement));
          kFunction(epsVoigt, bopI, bopJ, i, j, gp);
        }
      }
    }
  }

  template <typename ScalarType>
  auto calculateScalarImpl(const Requirement& par, ScalarAffordance affordance,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const -> ScalarType {
    if constexpr (FE::isMixed())
      return ScalarType{0.0};
    return energyFunction(par, dx)();
  }

  template <typename ScalarType>
  void calculateVectorImpl(const Requirement& par, VectorAffordance affordance,
                           typename Traits::template VectorType<ScalarType> force,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    if constexpr (FE::isMixed())
      return;
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    const auto eps          = strainFunction(par, dx);
    const auto fIntFunction = internalForcesFunction<ScalarType>(par, force, dx);

    // Internal forces
    for (const auto& [gpIndex, gp] : eps.viewOverIntegrationPoints()) {
      const auto epsVoigt = eps.evaluate(gpIndex, on(gridElement));
      const auto stresses = stress(epsVoigt);
      for (size_t i = 0; i < numberOfNodes_; ++i) {
        const auto bopI = eps.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        fIntFunction(stresses, bopI, i, gp);
      }
    }
  }
};

/**
 * \brief A helper function to create a linear elastic pre finite element.
 * \tparam MAT Type of the material.
 * \param mat Material parameters for the non-linear elastic element.
 * \return A linear elastic pre finite element.
 */
template <Concepts::GeometricallyLinearMaterial MAT>
auto linearElastic(const MAT& mat) {
  LinearElasticPre<MAT> pre(mat);

  return pre;
}
} // namespace Ikarus

#else
  #error LinearElastic depends on dune-localfefunctions, which is not included
#endif
