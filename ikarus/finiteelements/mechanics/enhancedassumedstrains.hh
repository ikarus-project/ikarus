// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file enhancedassumedstrains.hh
 * \brief Definition of the EAS class.
 * \ingroup mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/localfefunctions/derivativetransformators.hh>
  #include <dune/localfefunctions/linearAlgebraHelper.hh>
  #include <dune/localfefunctions/meta.hh>

  #include <ikarus/finiteelements/fehelper.hh>
  #include <ikarus/finiteelements/ferequirements.hh>
  #include <ikarus/finiteelements/mechanics/materials/tags.hh>
  #include <ikarus/finiteelements/mechanics/strainenhancements/easvariants.hh>
  #include <ikarus/utils/broadcaster/broadcastermessages.hh>
  #include <ikarus/utils/concepts.hh>

namespace Ikarus {

template <typename PreFE, typename FE, typename ES>
class EnhancedAssumedStrains;

/**
 * \brief A PreFE struct for Enhanced Assumed Strains.
 * \tparam ES The enhanced strain type.
 */
template <typename ES>
struct EnhancedAssumedStrainsPre
{
  int numberOfParameters{};

  template <typename PreFE, typename FE>
  using Skill = EnhancedAssumedStrains<PreFE, FE, ES>;
};

/**
 * \brief Wrapper class for using Enhanced Assumed Strains (EAS) with displacement based elements.
 *
 * \ingroup mechanics
 *
 * This class extends a displacement-based element to support Enhanced Assumed Strains.
 *
 * \tparam PreFE Type of the pre finite element.
 * \tparam FE Type of the finite element.
 * \tparam ESF The eas strain function.
 */
template <typename PreFE, typename FE, typename ESF>
class EnhancedAssumedStrains
    : public std::conditional_t<std::same_as<ESF, EAS::LinearStrain>,
                                ResultTypeBase<ResultTypes::linearStress, ResultTypes::linearStressFull>,
                                ResultTypeBase<ResultTypes::PK2Stress, ResultTypes::PK2StressFull>>
{
public:
  using Traits                 = PreFE::Traits;
  using Requirement            = FERequirements<FESolutions::displacement, FEParameter::loadfactor>;
  using LocalView              = typename Traits::LocalView;
  using Geometry               = typename Traits::Geometry;
  using GridView               = typename Traits::GridView;
  using Pre                    = EnhancedAssumedStrainsPre<ESF>;
  using EnhancedStrainFunction = ESF;

  static constexpr int myDim = Traits::mydim;

  template <typename ST>
  using VectorXOptRef = std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>;

  template <template <typename, int, int> class RT>
  using RTWrapperType = ResultWrapper<RT<typename Traits::ctype, Traits::mydim, Traits::worlddim>, ResultShape::Vector>;

  /**
   * \brief Constructor for Enhanced Assumed Strains elements.
   * \param pre The pre finite element
   */
  explicit EnhancedAssumedStrains(const Pre& pre) {
    static_assert(Concepts::Formulations::TotalLagrangian<FE::strainType, FE::stressType>,
                  "EAS method is only implemented for the total Lagrangian setting.");
    static_assert(
        not(FE::strainType == StrainTags::linear) or (std::same_as<EnhancedStrainFunction, EAS::LinearStrain>),
        "If FE::strainType is linear, then enhancedStrain must also be linear.");
    this->setEASType(pre.numberOfParameters);
  }

  /**
   * \brief Checks if the element is displacement-based and the EAS is turned off.
   *
   * \return True if the element is displacement-based, false otherwise.
   */
  bool isDisplacementBased() const { return easVariant_.isDisplacmentBased(); }

  /**
   * \brief Gets the variant representing the type of Enhanced Assumed Strains (EAS).
   *
   * \return Const reference to the EAS variant.
   */
  const auto& easVariant() const { return easVariant_; }

  /**
   * \brief Gets the number of EAS parameters based on the current EAS type.
   *
   * \return Number of EAS parameters.
   */
  auto numberOfInternalVariables() const { return easVariant_.numberOfInternalVariables(); }

  /**
   * \brief Calculates a requested result at a specific local position using the Enhanced Assumed Strains (EAS)
   * method.
   *
   * This function calculates the results at the specified local coordinates .
   * It takes into account the displacement-based element calculations and, if applicable, incorporates the EAS method
   * for enhanced accuracy.
   *
   * \param req The result requirements.
   * \param local The local coordinates at which results are to be calculated.
   * \return calculated result
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  requires(EnhancedAssumedStrains::template supportsResultType<RT>())
  auto calculateAtImpl(const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<2>) const {
    if constexpr (isSameResultType<RT, ResultTypes::linearStress> or isSameResultType<RT, ResultTypes::PK2Stress> or
                  isSameResultType<RT, ResultTypes::linearStressFull> or
                  isSameResultType<RT, ResultTypes::PK2StressFull>) {
      const auto ufunc     = underlying().displacementFunction(req);
      const auto rFunction = underlying().template resultFunction<RT>();
      auto disp            = Dune::viewAsFlatEigenVector(ufunc.coefficientsRef());
      const auto geo       = underlying().localView().element().geometry();

      RTWrapperType<RT> resultWrapper{};
      auto calculateAtContribution = [&]<typename EAST>(const EAST& easFunction) {
        Eigen::VectorXd alpha;
        alpha.setZero(numberOfInternalVariables());
        if constexpr (EAST::enhancedStrainSize != 0) {
          typename EAST::DType D;
          calculateDAndLMatrix(easFunction, req, D, L_);
          alpha = -D.inverse() * L_ * disp;
        }
        const auto enhancedStrain = EnhancedStrainFunction::value(geo, ufunc, local, easFunction, alpha);
        resultWrapper             = rFunction(enhancedStrain);
      };
      easVariant_(calculateAtContribution);
      return resultWrapper;
    }
    DUNE_THROW(Dune::NotImplemented, "The requested result type is not supported");
  }

  /**
   * \brief Sets the EAS type for 2D elements.
   *
   * \param numberOfInternalVariables The number of EAS parameters
   */
  void setEASType(int numberOfInternalVariables) {
    if (numberOfInternalVariables != 0)
      easApplicabilityCheck();
    easVariant_.setEASType(numberOfInternalVariables);
    initializeState();
  }

  /**
   * \brief Gets the internal state variable alpha for the EAS element.
   *
   * \return Internal state variable (alpha).
   */
  const auto& internalVariable() const { return alpha_; }

protected:
  void bindImpl() {
    assert(underlying().localView().bound());
    easVariant_.bind(underlying().localView().element().geometry());
    initializeState();
  }

  /**
   * \brief Updates the internal state variable alpha_ at the end of an iteration before the update of the displacements
   * done by the non-linear solver. See \cite bieberLockingHourglassingNonlinear2024 for implementation details and
   * further references.
   *
   * \param par The Requirement object.
   * \param correction The correction in displacement (DeltaD) vector passed based on which the internal state variable
   * alpha is to be updated.
   */
  void updateStateImpl(const Requirement& par,
                       const std::remove_reference_t<typename Traits::template VectorType<>>& correction) {
    using ScalarType = Traits::ctype;
    easApplicabilityCheck();
    auto correctAlpha = [&]<typename EAST>(const EAST& easFunction) {
      if constexpr (EAST::enhancedStrainSize != 0) {
        const auto& Rtilde      = calculateRtilde<ScalarType>(par);
        const auto localdxBlock = Ikarus::FEHelper::localSolutionBlockVector<Traits, Eigen::VectorXd, double>(
            correction, underlying().localView());
        const auto localdx               = Dune::viewAsFlatEigenVector(localdxBlock);
        constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
        Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize> D;
        calculateDAndLMatrix(easFunction, par, D, L_);
        this->alpha_ -= D.inverse() * (Rtilde + (L_ * localdx));
      }
    };

    easVariant_(correctAlpha);
  }

  inline void easApplicabilityCheck() const {
    const auto& numberOfNodes = underlying().numberOfNodes();
    if (not((numberOfNodes == 4 and Traits::mydim == 2) or (numberOfNodes == 9 and Traits::mydim == 2) or
            (numberOfNodes == 8 and Traits::mydim == 3)) and
        (not isDisplacementBased()))
      DUNE_THROW(Dune::NotImplemented, "EAS is only supported for Q1, Q2 and H1 elements");
  }

  template <typename ScalarType>
  void calculateMatrixImpl(const Requirement& par, const MatrixAffordance& affordance,
                           typename Traits::template MatrixType<> K,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    if (affordance != MatrixAffordance::stiffness)
      DUNE_THROW(Dune::NotImplemented, "MatrixAffordance not implemented: " + toString(affordance));
    easApplicabilityCheck();

    const auto geo           = underlying().localView().element().geometry();
    const auto& uFunction    = underlying().displacementFunction(par, dx);
    const auto& kMFunction   = underlying().template materialStiffnessMatrixFunction<ScalarType>(par, K, dx);
    const auto& kGFunction   = underlying().template geometricStiffnessMatrixFunction<ScalarType>(par, K, dx);
    const auto numberOfNodes = underlying().numberOfNodes();
    const auto& localBasis   = underlying().localBasis();

    auto calculateMatrixContribution = [&]<typename EAST>(const EAST& easFunction) {
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto enhancedStrain = EnhancedStrainFunction::value(geo, uFunction, gp.position(), easFunction, alpha_);
        const auto stresses       = underlying().stress(enhancedStrain);
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto E_dI = EnhancedStrainFunction::template firstDerivative<0>(geo, uFunction, localBasis, gpIndex,
                                                                                gp.position(), easFunction, alpha_, i);
          for (size_t j = 0; j < numberOfNodes; ++j) {
            const auto E_dJ = EnhancedStrainFunction::template firstDerivative<0>(
                geo, uFunction, localBasis, gpIndex, gp.position(), easFunction, alpha_, j);
            const auto E_dd = EnhancedStrainFunction::template secondDerivative<0>(
                geo, uFunction, localBasis, gpIndex, gp.position(), stresses, easFunction, alpha_, i, j);
            kMFunction(enhancedStrain, E_dI, E_dJ, i, j, gp);
            kGFunction(E_dd, i, j, gp);
          }
        }
      }

      if constexpr (EAST::enhancedStrainSize != 0) {
        typename EAST::DType D;
        calculateDAndLMatrix(easFunction, par, D, L_);
        K.template triangularView<Eigen::Upper>() -= L_.transpose() * D.inverse() * L_;
        K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
      }
    };
    easVariant_(calculateMatrixContribution);
  }

  template <typename ScalarType>
  inline ScalarType calculateScalarImpl(const Requirement& par, ScalarAffordance affordance,
                                        const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    easApplicabilityCheck();
    if (isDisplacementBased())
      return underlying().template energyFunction<ScalarType>(par, dx)();

    DUNE_THROW(Dune::NotImplemented,
               "EAS element do not support any scalar calculations, i.e. they are not derivable from a potential");
  }

  template <typename ScalarType>
  void calculateVectorImpl(const Requirement& par, VectorAffordance affordance,
                           typename Traits::template VectorType<ScalarType> force,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    if (affordance != VectorAffordance::forces)
      DUNE_THROW(Dune::NotImplemented, "VectorAffordance not implemented: " + toString(affordance));
    easApplicabilityCheck();

    const auto geo           = underlying().localView().element().geometry();
    const auto& uFunction    = underlying().displacementFunction(par, dx);
    const auto numberOfNodes = underlying().numberOfNodes();
    const auto& localBasis   = underlying().localBasis();
    const auto fIntFunction  = underlying().template internalForcesFunction<ScalarType>(par, force, dx);

    auto calculateForceContribution = [&]<typename EAST>(const EAST& easFunction) {
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto enhancedStrain = EnhancedStrainFunction::value(geo, uFunction, gp.position(), easFunction, alpha_);
        const auto stresses       = underlying().stress(enhancedStrain);
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto E_dI = EnhancedStrainFunction::template firstDerivative<0>(geo, uFunction, localBasis, gpIndex,
                                                                                gp.position(), easFunction, alpha_, i);
          fIntFunction(stresses, E_dI, i, gp);
        }
      }
      if constexpr (EAST::enhancedStrainSize != 0) {
        // Here, L_ and D doesn't see the ScalarType, while Rtilde sees it. This is needed because only then when we use
        // automatic differentiation w.r.t to the displacements, we get the correct static-condensed stiffness matrix.
        // If L_ and D sees the ScalarType, then their derivatives w.r.t. the displacements is also performed which is
        // then not equivalent to the static-condensed stiffness matrix in general.
        typename EAST::DType D;
        calculateDAndLMatrix(easFunction, par, D, L_);
        const auto& Rtilde = calculateRtilde(par, dx);
        force -= L_.transpose() * D.inverse() * Rtilde;
      }
    };
    easVariant_(calculateForceContribution);
  }

  template <typename MT, typename BC>
  void subscribeToImpl(BC& bc) {
    if constexpr (std::same_as<MT, NonLinearSolverMessages>) {
      using NLSState = typename BC::State;
      underlying().subscribe(bc, [&](NonLinearSolverMessages message, const NLSState& state) {
        if (message == NonLinearSolverMessages::CORRECTION_UPDATED)
          this->updateStateImpl(state.domain, state.correction);
      });
    }
  }

private:
  EAS::EASVariant<ESF, Geometry> easVariant_;
  mutable Eigen::MatrixXd L_;
  Eigen::VectorXd alpha_;

  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }

  /**
   * \brief Initializes the internal state variable alpha_ based on the number of EAS parameters.
   */
  void initializeState() {
    alpha_.resize(numberOfInternalVariables());
    alpha_.setZero();
  }

  template <int enhancedStrainSize>
  void calculateDAndLMatrix(const auto& easFunction, const auto& par,
                            Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>& DMat,
                            Eigen::MatrixX<double>& LMat) const {
    using namespace Dune;
    using namespace Dune::DerivativeDirections;

    const auto& uFunction    = underlying().displacementFunction(par);
    const auto geo           = underlying().localView().element().geometry();
    const auto numberOfNodes = underlying().numberOfNodes();
    const auto& localBasis   = underlying().localBasis();
    DMat.setZero();
    LMat.setZero(enhancedStrainSize, underlying().localView().size());
    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const auto enhancedStrain = EnhancedStrainFunction::value(geo, uFunction, gp.position(), easFunction, alpha_);
      const auto C              = underlying().materialTangent(enhancedStrain);
      const auto stresses       = underlying().stress(enhancedStrain);
      const double intElement   = geo.integrationElement(gp.position()) * gp.weight();
      const auto E_a  = EnhancedStrainFunction::template firstDerivative<1>(geo, uFunction, localBasis, gpIndex,
                                                                            gp.position(), easFunction, alpha_);
      const auto E_aa = EnhancedStrainFunction::template secondDerivative<1>(
          geo, uFunction, localBasis, gpIndex, gp.position(), stresses, easFunction, alpha_);
      DMat += (E_a.transpose() * C * E_a + E_aa) * intElement;
      for (size_t i = 0U; i < numberOfNodes; ++i) {
        const auto E_dI = EnhancedStrainFunction::template firstDerivative<0>(geo, uFunction, localBasis, gpIndex,
                                                                              gp.position(), easFunction, alpha_, i);
        const auto E_ad = EnhancedStrainFunction::template secondDerivative<2>(
            geo, uFunction, localBasis, gpIndex, gp.position(), stresses, easFunction, alpha_, i);
        LMat.template block<enhancedStrainSize, myDim>(0, myDim * i) +=
            (E_a.transpose() * C * E_dI + E_ad) * intElement;
      }
    }
  }

  template <typename ScalarType>
  Eigen::VectorX<ScalarType> calculateRtilde(const Requirement& par,
                                             const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    const auto geo         = underlying().localView().element().geometry();
    const auto& uFunction  = underlying().displacementFunction(par, dx);
    const auto& localBasis = underlying().localBasis();
    Eigen::VectorX<ScalarType> Rtilde;
    Rtilde.setZero(numberOfInternalVariables());

    auto calculateRtildeContribution = [&]<typename EAST>(const EAST& easFunction) {
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto enhancedStrain = EnhancedStrainFunction::value(geo, uFunction, gp.position(), easFunction, alpha_);
        const auto E_a = EnhancedStrainFunction::template firstDerivative<1>(geo, uFunction, localBasis, gpIndex,
                                                                             gp.position(), easFunction, alpha_);
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        auto stresses           = underlying().stress(enhancedStrain);
        Rtilde += (E_a.transpose() * stresses).eval() * intElement;
      }
    };

    easVariant_(calculateRtildeContribution);
    return Rtilde;
  }
};

/**
 * \brief A helper function to create an enhanced assumed strain pre finite element.
 * \tparam ES The strain tag that is enhanced.
 * \param numberOfInternalVariables Number of EAS parameters
 * \return An enhanced assumed strain pre finite element.
 */
template <typename ES = EAS::LinearStrain>
auto eas(int numberOfInternalVariables = 0) {
  EnhancedAssumedStrainsPre<ES> pre(numberOfInternalVariables);

  return pre;
}

} // namespace Ikarus

#else
  #error EnhancedAssumedStrains depends on dune-localfefunctions, which is not included
#endif
