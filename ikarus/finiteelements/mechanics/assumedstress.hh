// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file assumedstress.hh
 * \brief Definition of the AssumedStress class.
 * \ingroup mechanics
 */

#pragma once

#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/localfefunctions/derivativetransformators.hh>
  #include <dune/localfefunctions/linearAlgebraHelper.hh>
  #include <dune/localfefunctions/meta.hh>

  #include <ikarus/finiteelements/fehelper.hh>
  #include <ikarus/finiteelements/ferequirements.hh>
  #include <ikarus/finiteelements/mechanics/assumedstress/asfunctions.hh>
  #include <ikarus/finiteelements/mechanics/assumedstress/asvariants.hh>
  #include <ikarus/finiteelements/mechanics/materials/tags.hh>
  #include <ikarus/utils/broadcaster/broadcastermessages.hh>
  #include <ikarus/utils/concepts.hh>
  #include <ikarus/utils/linearalgebrahelper.hh>

namespace Ikarus {

template <typename PreFE, typename FE, typename ASType>
class AssumedStress;

/**
 * \brief A PreFE struct for Assumed Stress.
 * \tparam ES The stress type that is assumed.
 */
template <typename ASType>
struct AssumedStressPre
{
  int numberOfParameters{};

  template <typename PreFE, typename FE>
  using Skill = AssumedStress<PreFE, FE, ASType>;
};

/**
 * \brief Wrapper class for using Assumed Stress with displacement based elements.
 *
 * \ingroup mechanics
 *
 * This class extends a displacement-based element to support Assumed Stress.
 *
 * \tparam PreFE Type of the pre finite element.
 * \tparam FE Type of the finite element.
 * \tparam AST The as stress function.
 */
template <typename PreFE, typename FE, typename ASF>
class AssumedStress
    : public std::conditional_t<std::same_as<ASF, PS::LinearStress>,
                                ResultTypeBase<ResultTypes::linearStress, ResultTypes::linearStressFull>,
                                ResultTypeBase<ResultTypes::PK2Stress, ResultTypes::PK2StressFull>>

{
public:
  using Traits                = PreFE::Traits;
  using Requirement           = FERequirements<FESolutions::displacement, FEParameter::loadfactor>;
  using LocalView             = typename Traits::LocalView;
  using Geometry              = typename Traits::Geometry;
  using GridView              = typename Traits::GridView;
  using Pre                   = AssumedStressPre<ASF>;
  using AssumedStressFunction = ASF;

  template <typename ST>
  using VectorXOptRef = std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>;

  template <template <typename, int, int> class RT>
  using RTWrapperType = ResultWrapper<RT<typename Traits::ctype, Traits::mydim, Traits::worlddim>, ResultShape::Vector>;

  static constexpr int myDim     = Traits::mydim;
  static constexpr int strainDim = myDim * (myDim + 1) / 2;
  using StrainVector             = Eigen::Vector<double, strainDim>;
  using MaterialMatrix           = Eigen::Matrix<double, strainDim, strainDim>;

  /**
   * \brief Constructor for Assunmed Stress elements.
   * \param pre The pre finite element
   */
  explicit AssumedStress(const Pre& pre) {
    static_assert(Concepts::Formulations::TotalLagrangian<FE::strainType, FE::stressType>,
                  "PS method is only implemented for the total Lagrangian setting.");
    static_assert(not(FE::strainType == StrainTags::linear) or (std::same_as<ASF, PS::LinearStress>),
                  "If FE::strainType is linear, then the assumed stress must also be linear.");
    asVariant_.setAssumedStressType(pre.numberOfParameters);
  }

  /**
   * \brief Gets the variant representing the type of Assumed Stress.
   *
   * \return Const reference to the AssumedStress variant.
   */
  const auto& asVariant() const { return asVariant_; }

  /**
   * \brief Gets the number of AssumedStress parameters based on the current AssumedStress type.
   *
   * \return Number of AssumedStress parameters.
   */
  auto numberOfInternalVariables() const { return asVariant_.numberOfInternalVariables(); }

  /**
   * \brief Calculates a requested result at a specific local position using the Assumed Stress Method
   *
   * This function calculates the results at the specified local coordinates .
   * It takes into account the displacement-based element calculations and, incorporates the
   * AssumedStress method for enhanced accuracy.
   *
   * \param req The result requirements.
   * \param local The local coordinates at which results are to be calculated.
   * \return calculated result
   *
   * \tparam RT The type representing the requested result.
   */
  template <template <typename, int, int> class RT>
  requires(AssumedStress::template supportsResultType<RT>())
  auto calculateAtImpl(const Requirement& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<2>) const {
    if constexpr (isSameResultType<RT, ResultTypes::linearStress> or isSameResultType<RT, ResultTypes::PK2Stress> or
                  isSameResultType<RT, ResultTypes::linearStressFull> or
                  isSameResultType<RT, ResultTypes::PK2StressFull>) {
      const auto geo   = underlying().localView().element().geometry();
      const auto ufunc = underlying().displacementFunction(req);
      auto disp        = Dune::viewAsFlatEigenVector(ufunc.coefficientsRef());

      RTWrapperType<RT> resultWrapper{};
      auto calculateAtContribution = [&]<typename AssumedStressT>(const AssumedStressT& asFunction) {
        Eigen::VectorXd beta;
        beta.setZero(numberOfInternalVariables());
        const auto& Rtilde = calculateRtilde<double>(req);
        typename AssumedStressT::HType H;
        calculateHAndGMatrix(asFunction, req, H, G_);
        beta = H.inverse() * (G_ * disp);

        const auto SVoigt = AssumedStressFunction::value(geo, local, asFunction, beta);

        if constexpr ((isSameResultType<RT, ResultTypes::PK2StressFull> or
                       isSameResultType<RT, ResultTypes::linearStressFull>) and
                      requires { underlying().material().underlying(); }) {
          resultWrapper = enlargeStressAnsatz(SVoigt);
        } else {
          resultWrapper = SVoigt;
        }
      };
      asVariant_(calculateAtContribution);
      return resultWrapper;
    }
    DUNE_THROW(Dune::NotImplemented, "The requested result type is not supported");
  }

  /**
   * \brief Sets the AssumedStress type for 2D elements.
   *
   * \param numberOfInternalVariables The number of AssumedStress parameters
   */
  void setAssumedStressType(int numberOfInternalVariables) {
    asApplicabilityCheck();
    asVariant_.setAssumedStressType(numberOfInternalVariables);
    initializeState();
  }

  /**
   * \brief Gets the internal state variable beta for the AssumedStress element.
   *
   * \return Internal state variable (beta).
   */
  const auto& internalVariable() const { return beta_; }

protected:
  void bindImpl() {
    assert(underlying().localView().bound());
    asVariant_.bind(underlying().localView().element().geometry());
    initializeState();
  }

  /**
   * \brief Updates the internal state variable beta_ at the end of an iteration before the update of the displacements
   * done by the non-linear solver. See \cite bieberLockingHourglassingNonlinear2024 for implementation details and
   * further references.
   *
   * \param par The Requirement object.
   * \param correction The correction in displacement (DeltaD) vector passed based on which the internal state variable
   * beta is to be updated.
   */
  void updateStateImpl(const Requirement& par,
                       const std::remove_reference_t<typename Traits::template VectorType<>>& correction) {
    using ScalarType = Traits::ctype;
    asApplicabilityCheck();
    auto correctbeta = [&]<typename AssumedStressT>(const AssumedStressT& asFunction) {
      if (correction.size() != par.globalSolution().size())
        DUNE_THROW(Dune::NotImplemented,
                   "Solution vector and correction vector should be of the same size. Check if DBCOption::Full is "
                   "used. The sizes are " +
                       std::to_string(par.globalSolution().size()) + " and " + std::to_string(correction.size()));
      const auto localdxBlock = Ikarus::FEHelper::localSolutionBlockVector<Traits, Eigen::VectorXd, double>(
          correction, underlying().localView());
      const auto localdx = Dune::viewAsFlatEigenVector(localdxBlock);

      typename AssumedStressT::HType H;
      calculateHAndGMatrix(asFunction, par, H, G_);
      const auto& Rtilde = calculateRtilde<ScalarType>(par);
      this->beta_ += H.inverse() * (Rtilde + (G_ * localdx));
    };

    asVariant_(correctbeta);
    calculateMaterialInversion();
  }

  inline void asApplicabilityCheck() const {
    const auto& numberOfNodes = underlying().numberOfNodes();
    if (not((numberOfNodes == 4 and Traits::mydim == 2) or (numberOfNodes == 8 and Traits::mydim == 3)))
      DUNE_THROW(Dune::NotImplemented, "AssumedStress is only supported for Q1 and H1 elements" +
                                           std::to_string(numberOfNodes) + std::to_string(Traits::mydim));
  }

  template <typename ScalarType>
  void calculateMatrixImpl(const Requirement& par, const MatrixAffordance& affordance,
                           typename Traits::template MatrixType<> K,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    using namespace Dune;
    using namespace Dune::DerivativeDirections;

    if (affordance != MatrixAffordance::stiffness)
      DUNE_THROW(Dune::NotImplemented, "MatrixAffordance not implemented: " + toString(affordance));
    asApplicabilityCheck();

    const auto geo             = underlying().localView().element().geometry();
    const auto& strainFunction = underlying().strainFunction(par, dx);
    const auto& kGFunction     = underlying().template geometricStiffnessMatrixFunction<ScalarType>(par, K, dx);
    const auto numberOfNodes   = underlying().numberOfNodes();
    const auto& localBasis     = underlying().localBasis();

    auto calculateMatrixContribution = [&]<typename AssumedStressT>(const AssumedStressT& asFunction) {
      for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
        const auto SVoigt = AssumedStressFunction::value(geo, gp.position(), asFunction, beta_);

        for (size_t i = 0; i < numberOfNodes; ++i)
          for (size_t j = 0; j < numberOfNodes; ++j) {
            const auto E_dd =
                strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i, j)), along(SVoigt), on(gridElement));
            kGFunction(E_dd, i, j, gp);
          }
      }

      typename AssumedStressT::HType H;
      calculateHAndGMatrix(asFunction, par, H, G_, dx);

      K.template triangularView<Eigen::Upper>() += G_.transpose() * H.inverse() * G_;
      K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
    };
    asVariant_(calculateMatrixContribution);
  }

  template <typename ScalarType>
  inline ScalarType calculateScalarImpl(const Requirement& par, ScalarAffordance affordance,
                                        const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    asApplicabilityCheck();
    ScalarType energy = 0.0;

    // This will not work in the context of autodiff because `AutoDiffFE` class doesn't include static condensation of
    // the internal variable
    if constexpr (not Concepts::AutodiffScalar<ScalarType>) {
      const auto geo            = underlying().localView().element().geometry();
      const auto strainFunction = underlying().strainFunction(par, dx);

      auto calculateScalarContribution = [&]<typename AssumedStressT>(const AssumedStressT& asFunction) {
        for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
          const auto SVoigt = AssumedStressFunction::value(geo, gp.position(), asFunction, beta_);
          const auto& Es    = strainStates_[gpIndex];
          const auto EVoigt = strainFunction.evaluate(gp.position(), on(gridElement));

          energy += ((SVoigt.transpose() * EVoigt - ScalarType(0.5) * (SVoigt.transpose() * Es)).eval() *
                     geo.integrationElement(gp.position()) * gp.weight())
                        .eval()[0];
        }
      };

      asVariant_(calculateScalarContribution);
    } else {
      DUNE_THROW(Dune::NotImplemented,
                 "AssumedStress element does not support scalar calculations for autodiff scalars");
    }
    return energy;
  }

  template <typename ScalarType>
  void calculateVectorImpl(const Requirement& par, VectorAffordance affordance,
                           typename Traits::template VectorType<ScalarType> force,
                           const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;

    if (affordance != VectorAffordance::forces)
      DUNE_THROW(Dune::NotImplemented, "VectorAffordance not implemented: " + toString(affordance));
    asApplicabilityCheck();

    const auto geo             = underlying().localView().element().geometry();
    const auto& uFunction      = underlying().displacementFunction(par, dx);
    const auto& strainFunction = underlying().strainFunction(par, dx);
    const auto numberOfNodes   = underlying().numberOfNodes();
    const auto& localBasis     = underlying().localBasis();
    const auto fIntFunction    = underlying().template internalForcesFunction<ScalarType>(par, force, dx);

    auto calculateForceContribution = [&]<typename AssumedStressT>(const AssumedStressT& asFunction) {
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto SVoigt = AssumedStressFunction::value(geo, gp.position(), asFunction, beta_);

        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto E_dI = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          fIntFunction(SVoigt, E_dI, i, gp);
        }
      }
      typename AssumedStressT::HType H;
      calculateHAndGMatrix(asFunction, par, H, G_);
      const auto Rtilde = calculateRtilde(par, dx);
      force += G_.transpose() * H.inverse() * Rtilde;
    };
    asVariant_(calculateForceContribution);
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
  PS::AssumedStressVariant<AssumedStressFunction, Geometry> asVariant_;
  mutable Eigen::MatrixXd G_;
  Eigen::VectorXd beta_;
  std::vector<StrainVector> strainStates_;
  std::vector<MaterialMatrix> invertedMaterialStates_;

  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }

  /**
   * \brief Initializes the internal state variable beta_ based on the number of AssumedStress parameters.
   */
  void initializeState() {
    beta_.resize(numberOfInternalVariables());
    beta_.setZero();
    strainStates_.resize(underlying().localBasis().integrationPointSize(), StrainVector::Zero());
    invertedMaterialStates_.resize(underlying().localBasis().integrationPointSize(), MaterialMatrix::Zero());
    calculateMaterialInversion();
  }

  template <typename ScalarType, int assumedStressSize>
  void calculateHAndGMatrix(const auto& asFunction, const auto& par,
                            Eigen::Matrix<double, assumedStressSize, assumedStressSize>& H,
                            Eigen::MatrixX<ScalarType>& G, const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;

    const auto& uFunction      = underlying().displacementFunction(par);
    const auto& strainFunction = underlying().strainFunction(par, dx);
    const auto geo             = underlying().localView().element().geometry();
    const auto numberOfNodes   = underlying().numberOfNodes();
    const auto& localBasis     = underlying().localBasis();
    H.setZero();
    G.setZero(assumedStressSize, underlying().localView().size());
    for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
      const double intElement = geo.integrationElement(gp.position()) * gp.weight();

      const auto SVoigt = AssumedStressFunction::value(geo, gp.position(), asFunction, beta_);
      const auto& D     = invertedMaterialStates_[gpIndex];

      const auto S_b  = AssumedStressFunction::template firstDerivative(geo, uFunction, localBasis, gpIndex,
                                                                        gp.position(), asFunction, beta_);
      const auto S_bb = AssumedStressFunction::template secondDerivative(geo, uFunction, localBasis, gpIndex,
                                                                         gp.position(), SVoigt, asFunction, beta_);

      H += (S_b.transpose() * D * S_b + S_bb) * intElement;

      for (size_t i = 0U; i < numberOfNodes; ++i) {
        const auto E_dI = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));

        G.template block<assumedStressSize, Traits::worlddim>(0, myDim * i) += S_b.transpose() * E_dI * intElement;
      }
    }
  }

  template <typename ScalarType>
  Eigen::VectorX<ScalarType> calculateRtilde(const Requirement& par,
                                             const VectorXOptRef<ScalarType>& dx = std::nullopt) const {
    using namespace Dune;
    using namespace Dune::DerivativeDirections;

    const auto geo            = underlying().localView().element().geometry();
    const auto& uFunction     = underlying().displacementFunction(par, dx);
    const auto strainFunction = underlying().strainFunction(par, dx);
    const auto& localBasis    = underlying().localBasis();

    Eigen::VectorX<ScalarType> Rtilde;
    Rtilde.setZero(numberOfInternalVariables());

    auto calculateRtildeContribution = [&]<typename AssumedStressT>(const AssumedStressT& asFunction) {
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto SVoigt = AssumedStressFunction::value(geo, gp.position(), asFunction, beta_);
        const auto S_b    = AssumedStressFunction::template firstDerivative(geo, uFunction, localBasis, gpIndex,
                                                                            gp.position(), asFunction, beta_);

        const auto& Es    = strainStates_[gpIndex];
        const auto EVoigt = strainFunction.evaluate(gp.position(), on(gridElement));

        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        Rtilde += (S_b.transpose() * (EVoigt - Es)).eval() * intElement;
      }
    };

    asVariant_(calculateRtildeContribution);
    return Rtilde;
  }

  void calculateMaterialInversion() {
    auto saveMaterialState = [&]<typename AssumedStressT>(const AssumedStressT& asFunction) {
      const auto geo = underlying().localView().element().geometry();

      for (const auto& [gpIndex, gp] : underlying().localBasis().viewOverIntegrationPoints()) {
        const auto SVoigt = AssumedStressFunction::value(geo, gp.position(), asFunction, beta_);
        const auto& EsOld = strainStates_[gpIndex];
        std::tie(invertedMaterialStates_[gpIndex], strainStates_[gpIndex]) =
            underlying().material().template materialInversion<FE::strainType, true>(SVoigt, EsOld);
      }
    };
    asVariant_(saveMaterialState);
  }

  auto enlargeStressAnsatz(const StrainVector& SVoigt, size_t qpIndexForIteration = 0) const
  requires(decltype(underlying().material())::isReduced)
  {
    auto [_, Es] = underlying().material().template materialInversion<FE::strainType, true>(
        SVoigt, strainStates_[qpIndexForIteration]);
    auto Esfull = enlargeIfReduced<decltype(underlying().material())>(Es).eval();
    return underlying().material().underlying().template stresses<FE::strainType, true>(Esfull);
  }
};

/**
 * \brief A helper function to create an assumed stress pre finite element.
 * \param numberOfInternalVariables Number of AssumedStress parameters
 * \return An assumed stress pre finite element.
 */
template <typename ASType = PS::LinearStress>
auto assumedStress(int numberOfInternalVariables) {
  return AssumedStressPre<ASType>(numberOfInternalVariables);
}

} // namespace Ikarus

#else
  #error AssumedStress depends on dune-localfefunctions, which is not included
#endif
