// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file enhancedassumedstrains.hh
 * \brief Definition of the EAS class.
 * \ingroup  mechanics
 */

#pragma once
#if HAVE_DUNE_LOCALFEFUNCTIONS
  #include <dune/localfefunctions/derivativetransformators.hh>
  #include <dune/localfefunctions/linearAlgebraHelper.hh>
  #include <dune/localfefunctions/meta.hh>

  #include <ikarus/finiteelements/fehelper.hh>
  #include <ikarus/finiteelements/ferequirements.hh>
  #include <ikarus/finiteelements/mechanics/easvariants.hh>
  #include <ikarus/utils/eigendunetransformations.hh>

namespace Ikarus {

/**
 * \brief Wrapper class for using Enhanced Assumed Strains (EAS) with displacement based elements.
 *
 * \ingroup mechanics
 *
 * This class extends a displacement-based element to support Enhanced Assumed Strains.
 *
 * \tparam DFE The base displacement-based element type.
 */
template <typename DFE>
class EnhancedAssumedStrains : public DFE
{
public:
  using DisplacementBasedElement = DFE;
  using FERequirementType        = typename DisplacementBasedElement::FERequirementType;
  using LocalView                = typename DisplacementBasedElement::LocalView;
  using Geometry                 = typename DisplacementBasedElement::Geometry;
  using GridView                 = typename DisplacementBasedElement::GridView;
  using Traits                   = typename DisplacementBasedElement::Traits;

  using DisplacementBasedElement::localView;

  /**
* \brief Constructor for Enhanced Assumed Strains elements.
  * \details Disabling this forwarding constructor if the argument provided is EnhancedAssumedStrains itself, to
forward the
  // calls to the implicit copy constructor
* \tparam Args Variadic template for constructor arguments.
* \param args Constructor arguments forwarded to the base class.
*/
  template <typename... Args>
  requires(
      not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, EnhancedAssumedStrains>)
  explicit EnhancedAssumedStrains(Args&&... args)
      : DisplacementBasedElement(std::forward<Args>(args)...) {}

  /**
   * \brief Calculates a scalar quantity for the element.
   *
   * This function calculates a scalar quantity for the element based on the FERequirementType.
   *
   * \param par The FERequirementType object.
   * \return Computed scalar quantity.
   */
  double calculateScalar(const FERequirementType& par) const {
    if (isDisplacementBased())
      return DisplacementBasedElement::calculateScalar(par);
    DUNE_THROW(Dune::NotImplemented,
               "EAS element do not support any scalar calculations, i.e. they are not derivable from a potential");
  }

  /**
   * \brief Checks if the element is displacement-based and the EAS is turned off.
   *
   * \return True if the element is displacement-based, false otherwise.
   */
  bool isDisplacementBased() const { return easVariant_.isDisplacmentBased(); }

  /**
   * \brief Calculates vectorial quantities for the element.
   *
   * This function calculates the vectorial quantities for the element based on the FERequirementType.
   *
   * \param par The FERequirementType object.
   * \param force Vector to store the calculated forces.
   */
  inline void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> force) const {
    calculateVectorImpl<double>(par, force);
  }

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
  auto numberOfEASParameters() const { return easVariant_.numberOfEASParameters(); }

  /**
   * \brief Calculates the matrix for the element.
   *
   * This function calculates the matrix for the element based on the FERequirementType.
   *
   * \param par The FERequirementType object.
   * \param K Matrix to store the calculated stiffness.
   */
  void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> K) const {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;

    /// fill K with displacement-based stiffness.
    /// It is assumed to be assembled block-wise on element level.
    /// This means the displacements x,y,z of node I are grouped together.
    DisplacementBasedElement::calculateMatrix(par, K);

    auto calculateMatrixContribution = [&]<typename EAST>(const EAST& easFunction) {
      typename EAST::DType D;
      calculateDAndLMatrix(easFunction, par, D, L_);

      K.template triangularView<Eigen::Upper>() -= L_.transpose() * D.inverse() * L_;
      K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
    };
    easVariant_(calculateMatrixContribution);
  }

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
  requires(DisplacementBasedElement::template canProvideResultType<RT>())
  auto calculateAt(const FERequirementType& req, const Dune::FieldVector<double, Traits::mydim>& local) const {
    if (isDisplacementBased())
      return DisplacementBasedElement::template calculateAt<RT>(req, local);

    using RTWrapper = ResultWrapper<RT<typename Traits::ctype, Traits::mydim, Traits::worlddim>, ResultShape::Vector>;
    using namespace Dune::Indices;
    using namespace Dune::DerivativeDirections;

    if constexpr (isSameResultType<RT, ResultType::linearStress>) {
      const auto C              = DisplacementBasedElement::materialTangentFunction(req);
      const auto& numberOfNodes = DisplacementBasedElement::numberOfNodes();
      auto uFunction            = DisplacementBasedElement::displacementFunction(req);
      const auto disp           = Dune::viewAsFlatEigenVector(uFunction.coefficientsRef());

      return std::visit(
          [&]<typename EAST>(const EAST& easFunction) {
            if constexpr (not std::is_same_v<std::monostate, EAST>) {
              constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
              Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize> D;
              calculateDAndLMatrix(easFunction, req, D, L_);
              const auto alpha = (-D.inverse() * L_ * disp).eval();
              const auto M     = easFunction.calcM(local);
              const auto CEval = C(local);
              auto easStress   = (CEval * M * alpha).eval();
              return RTWrapper{DisplacementBasedElement::template calculateAt<RT>(req, local).asVec() + easStress};
            } else
              return DisplacementBasedElement::template calculateAt<RT>(req, local);
          },
          easVariant_);
    } else
      static_assert(Dune::AlwaysFalse<DFE>::value, "The requested result type is NOT implemented.");
    return result;
    const auto C              = DisplacementBasedElement::materialTangentFunction(req);
    const auto& numberOfNodes = DisplacementBasedElement::numberOfNodes();
    auto uFunction            = DisplacementBasedElement::displacementFunction(req);

    auto calculateAtContribution = [&]<typename EAST>(const EAST& easFunction) {
      typename EAST::DType D;
      calculateDAndLMatrix(easFunction, req, D, L_);

      const auto disp  = Dune::viewAsFlatEigenVector(uFunction.coefficientsRef());
      const auto alpha = (-D.inverse() * L_ * disp).eval();
      const auto M     = easFunction.calcM(local);
      const auto CEval = C(local);
      auto easStress   = (CEval * M * alpha).eval();

      resultVector.resize(3, 1);
      resultVector = resultVector + easStress;
    };
    easVariant_(calculateAtContribution);
    return resultVector;
  }

  /**
   * \brief Sets the EAS type for 2D elements.
   *
   * \param numberOfEASParameters The number of EAS parameters
   */
  void setEASType(int numberOfEASParameters) {
    const auto& numberOfNodes = DisplacementBasedElement::numberOfNodes();
    if (not((numberOfNodes == 4 and Traits::mydim == 2) or (numberOfNodes == 8 and Traits::mydim == 3)) and
        (not isDisplacementBased()))
      DUNE_THROW(Dune::NotImplemented, "EAS only supported for Q1 or H1 elements");
    easVariant_.setEASType(numberOfEASParameters, localView().element().geometry());
  }

protected:
  template <typename ScalarType>
  inline ScalarType calculateScalarImpl(
      const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
    if (isDisplacementBased())
      return DisplacementBasedElement::template calculateScalarImpl<ScalarType>(par, dx);
    DUNE_THROW(Dune::NotImplemented,
               "EAS element do not support any scalar calculations, i.e. they are not derivable from a potential");
  }

  template <typename ScalarType>
  void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                           const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
    using namespace Dune;
    using namespace Dune::DerivativeDirections;

    DisplacementBasedElement::calculateVectorImpl(par, force, dx);

    const auto uFunction      = DisplacementBasedElement::displacementFunction(par, dx);
    auto strainFunction       = DisplacementBasedElement::strainFunction(par, dx);
    const auto& numberOfNodes = DisplacementBasedElement::numberOfNodes();

    auto C = DisplacementBasedElement::materialTangentFunction(par);

    auto calculateForceContribution = [&]<typename EAST>(const EAST& easFunction) {
      typename EAST::DType D;
      calculateDAndLMatrix(easFunction, par, D, L_);

      const auto disp  = Dune::viewAsFlatEigenVector(uFunction.coefficientsRef());
      const auto alpha = (-D.inverse() * L_ * disp).eval();
      const auto geo   = localView().element().geometry();

      for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
        const auto M            = easFunction.calcM(gp.position());
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        const auto CEval        = C(gpIndex);
        auto stresses           = (CEval * M * alpha).eval();
        for (size_t i = 0; i < numberOfNodes; ++i) {
          const auto bopI = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          force.template segment<Traits::worlddim>(Traits::worlddim * i) += bopI.transpose() * stresses * intElement;
        }
      }
    };
    easVariant_(calculateForceContribution);
  }

private:
  EAS::EASVariant<Geometry> easVariant_;
  mutable Eigen::MatrixXd L_;

  template <int enhancedStrainSize>
  void calculateDAndLMatrix(const auto& easFunction, const auto& par,
                            Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>& DMat,
                            Eigen::MatrixXd& LMat) const {
    using namespace Dune;
    using namespace Dune::DerivativeDirections;

    auto strainFunction      = DisplacementBasedElement::strainFunction(par);
    const auto C             = DisplacementBasedElement::materialTangentFunction(par);
    const auto geo           = localView().element().geometry();
    const auto numberOfNodes = DisplacementBasedElement::numberOfNodes();
    DMat.setZero();
    LMat.setZero(enhancedStrainSize, localView().size());
    for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
      const auto M      = easFunction.calcM(gp.position());
      const auto CEval  = C(gpIndex);
      const double detJ = geo.integrationElement(gp.position());
      DMat += M.transpose() * CEval * M * detJ * gp.weight();
      for (size_t i = 0U; i < numberOfNodes; ++i) {
        const size_t I = Traits::worlddim * i;
        const auto Bi  = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
        LMat.template block<enhancedStrainSize, Traits::worlddim>(0, I) +=
            M.transpose() * CEval * Bi * detJ * gp.weight();
      }
    }
  }
};
} // namespace Ikarus

#else
  #error EnhancedAssumedStrains depends on dune-localfefunctions, which is not included
#endif
