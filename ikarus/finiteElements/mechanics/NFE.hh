// SPDX-FileCopyrightText: 2022 The dune-localfefunction developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <dune/localfefunctions/impl/clonableLocalFunction.hh>

#include <concepts>

#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/localFunctionHelper.hh>
#include <dune/localfefunctions/localFunctionInterface.hh>

#include <Eigen/Core>
#include <Eigen/Dense>
//#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Dune {

  template <typename DuneBasis_, typename CoeffContainer_, typename Geometry_, std::size_t ID = 0,
            typename LinAlg = Dune::DefaultLinearAlgebra>
  class EmbeddedLocalFunction
      : public LocalFunctionInterface<EmbeddedLocalFunction<DuneBasis_, CoeffContainer_, Geometry_, ID, LinAlg>>,
        public ClonableLocalFunction<EmbeddedLocalFunction<DuneBasis_, CoeffContainer_, Geometry_, ID, LinAlg>> {
    using Interface = LocalFunctionInterface<EmbeddedLocalFunction>;

    template <size_t ID_ = 0>
    static constexpr int orderID = ID_ == ID ? 1 : 0;

  public:
    using DuneBasis =DuneBasis_;
    using CoeffContainer =CoeffContainer_;
    using Geometry =Geometry_;
    friend Interface;
    friend ClonableLocalFunction<EmbeddedLocalFunction>;
    constexpr EmbeddedLocalFunction(
        const Dune::CachedLocalBasis<DuneBasis>& p_basis, const CoeffContainer& coeffs_,
        const std::shared_ptr<const Geometry>& geo,
        Dune::template index_constant<ID> = Dune::template index_constant<std::size_t(0)>{})
        : basis_{p_basis},
          coeffs{coeffs_},
          geometry_{geo}  //          ,coeffsAsMat{Dune::viewAsEigenMatrixFixedDyn(coeffs)}
    {}

    using Traits = LocalFunctionTraits<EmbeddedLocalFunction>;

    static constexpr bool isLeaf = true;
    static constexpr std::array<int, 1> id{ID};

    using LinearAlgebra                                     = LinAlg;
    static constexpr bool providesDerivativeTransformations = true;

    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateDerivativeImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                       const LocalFunctionEvaluationArgs_& localFunctionArgs);

    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateFunctionImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                     const LocalFunctionEvaluationArgs_& localFunctionArgs);

    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    //    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = Traits::valueSize;
    //    /** \brief Dimension of the coeffs */
    static constexpr int correctionSize = Traits::correctionSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Traits::gridDim;
    /** \brief Dimension of the world where this function is mapped to from the reference element */
    static constexpr int worldDimension = Traits::worldDimension;
    /** \brief Type for coordinate vector in world space */
    using FunctionReturnType = typename Traits::FunctionReturnType;
    /** \brief Type for the directional derivatives */
    using AlongType = typename Traits::AlongType;
    /** \brief The manifold where the function values lives in */
    using Manifold = typename Traits::Manifold;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = typename Traits::Jacobian;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename Traits::JacobianColType;
    /** \brief Type for the derivatives wrT the coefficients in the local tangent space base*/
    using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
    /** \brief Type for the derivatives wrT the coefficients in the embedding space*/
    using CoeffDerivEukMatrix = typename Traits::CoeffDerivEukMatrix;
    /** \brief Type for the derivatives wrT the coefficients in the embedding space on the left and in the tangent space
     * basis on the right*/
    using CoeffDerivEukRieMatrix =
        typename DefaultLinearAlgebra::template FixedSizedMatrix<ctype, valueSize, correctionSize>;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

    const auto& coefficientsRef() const { return coeffs; }
    auto& coefficientsRef() { return coeffs; }
    auto& geometry() const { return geometry_; }

    const Dune::CachedLocalBasis<DuneBasis>& basis() const { return basis_; }

    template <typename OtherType>
    struct rebind {
      using other = EmbeddedLocalFunction<
          DuneBasis, typename Std::Rebind<CoeffContainer, typename Manifold::template rebind<OtherType>::other>::other,
          Geometry, ID, LinAlg>;
    };

  private:
    static auto tryToCallDerivativeOfProjectionWRTposition(const FunctionReturnType& valE) {
      if constexpr (requires { Manifold::derivativeOfProjectionWRTposition(valE); })
        return Manifold::derivativeOfProjectionWRTposition(valE);
      else
        static_assert(
            requires { Manifold::derivativeOfProjectionWRTposition(valE); },
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    static auto tryToCallSecondDerivativeOfProjectionWRTposition(const FunctionReturnType& valE,
                                                                 const AlongType& along) {
      if constexpr (requires { Manifold::secondDerivativeOfProjectionWRTposition(valE, along); })
        return Manifold::secondDerivativeOfProjectionWRTposition(valE, along);
      else
        static_assert(
            requires { Manifold::secondDerivativeOfProjectionWRTposition(valE, along); },
            " Your passed manifold does not implement secondDerivativeOfProjectionWRTposition.");
    }

    static auto tryToCallThirdDerivativeOfProjectionWRTposition(const FunctionReturnType& valE, const AlongType& along,
                                                                const AlongType& along2) {
      if constexpr (requires { Manifold::thirdDerivativeOfProjectionWRTposition(valE, along, along2); })
        return Manifold::thirdDerivativeOfProjectionWRTposition(valE, along, along2);
      else
        static_assert(
            requires { Manifold::thirdDerivativeOfProjectionWRTposition(valE, along, along2); },
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const On<TransformArgs...>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
//      FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      return  J;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                         int spaceIndex, const On<TransformArgs...>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      JacobianColType Jcol    = evaluateEmbeddingJacobianColImpl(dNTransformed, spaceIndex);
      return  Jcol;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           int coeffsIndex, const On<TransformArgs...>&) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      return N[coeffsIndex] * coeffs[coeffsIndex].orthonormalFrame();
    }


    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           const std::array<size_t, 2>& coeffsIndex,
                                                           const Along<AlongArgs...>& alongArgs,
                                                           const On<TransformArgs...>&) const {
      const auto& N                 = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      CoeffDerivEukMatrix ddt;
      ddt.setZero();

      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const CoeffDerivEukMatrix dt = N[coeffsIndex[0]]*CoeffDerivEukMatrix::Identity();
        auto& unitVec                = coeffs[coeffsIndex[0]];
        auto& unitVecVal             = unitVec.getValue();
        auto scal                    = inner(unitVecVal, dt * std::get<0>(alongArgs.args));
        auto idmat                   = createScaledIdentityMatrix<ctype, valueSize, valueSize>(scal);
        ddt -= idmat;
      }

      return transposeEvaluated(coeffs[coeffsIndex[0]].orthonormalFrame()) * ddt
             * coeffs[coeffsIndex[1]].orthonormalFrame();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    std::array<CoeffDerivEukRieMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const On<TransformArgs...>& transArgs) const {
      std::array<CoeffDerivEukRieMatrix, gridDim> WarrayRie;
      const auto BLA = coeffs[coeffsIndex].orthonormalFrame();
      for (int dir = 0; dir < gridDim; ++dir)
        WarrayRie[dir] = coeff(dNTransformed, coeffsIndex, dir) * BLA;

      return WarrayRie;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
        const On<TransformArgs...>& transArgs) const {
      return coeff(dNTransformed, coeffsIndex, spatialIndex) * coeffs[coeffsIndex].orthonormalFrame();
    }


    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const Along<AlongArgs...>& alongArgs, const On<TransformArgs...>& transArgs) const {
      CoeffDerivEukMatrix ChiArrayEuk;
      setZero(ChiArrayEuk);

      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const auto& along             = std::get<0>(alongArgs.args);
        for (int i = 0; i < gridDim; ++i) {
          const CoeffDerivEukMatrix W =coeff(dNTransformed, coeffsIndex[0], i)*CoeffDerivEukMatrix::Identity();
          ChiArrayEuk -= createScaledIdentityMatrix<ctype, valueSize, valueSize>(
              inner(coeffs[coeffsIndex[0]].getValue(), W * col(along, i)));
        }
      }
      CoeffDerivMatrix ChiArrayRie;
      const auto BLA0T = eval(transposeEvaluated(coeffs[coeffsIndex[0]].orthonormalFrame()));
      const auto BLA1  = coeffs[coeffsIndex[1]].orthonormalFrame();
      ChiArrayRie      = BLA0T * ChiArrayEuk * BLA1;
      return ChiArrayRie;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const int spatialIndex, const Along<AlongArgs...>& alongArgs, const On<TransformArgs...>& transArgs) const {
//      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
//      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
//      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
//      const Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      const auto& along             = std::get<0>(alongArgs.args);
//      const CoeffDerivEukMatrix S   = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along);
//      const auto chi          = tryToCallThirdDerivativeOfProjectionWRTposition(valE, along, col(J, spatialIndex));
//      const auto& NI          = N[coeffsIndex[0]];
//      const auto& NJ          = N[coeffsIndex[1]];
//      const auto& dNIdi       = coeff(dNTransformed, coeffsIndex[0], spatialIndex);
//      const auto& dNJdi       = coeff(dNTransformed, coeffsIndex[1], spatialIndex);
      CoeffDerivEukMatrix Chi;
      Chi.setZero();

      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const CoeffDerivEukMatrix W = coeff(dNTransformed, coeffsIndex[0], spatialIndex)*CoeffDerivEukMatrix::Identity();
        Chi -= createScaledIdentityMatrix<ctype, valueSize, valueSize>(
            inner(coeffs[coeffsIndex[0]].getValue(), W * along));
      }
      return (transposeEvaluated(coeffs[coeffsIndex[0]].orthonormalFrame()) * Chi
              * coeffs[coeffsIndex[1]].orthonormalFrame());
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    FunctionReturnType evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                            [[maybe_unused]] const On<TransformArgs...>&) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      return evaluateEmbeddingFunctionImpl(N);
    }

    JacobianColType evaluateEmbeddingJacobianColImpl(const AnsatzFunctionJacobian& dN, int spaceIndex) const {
      JacobianColType Jcol;
      setZero(Jcol);
      for (size_t j = 0; j < Rows<JacobianColType>::value; ++j) {
        for (size_t i = 0; i < coeffs.size(); ++i)
          Jcol[j] += coeffs[i].getValue()[j] * coeff(dN, i, spaceIndex);
      }

      return Jcol;
    }

    Jacobian evaluateEmbeddingJacobianImpl(const AnsatzFunctionJacobian& dN) const {
      Jacobian J;
      setZero(J);
      for (size_t j = 0; j < gridDim; ++j)
        for (size_t k = 0; k < valueSize; ++k)
          for (size_t i = 0; i < coeffs.size(); ++i)
            coeff(J, k, j) += coeffs[i].getValue()[k] * coeff(dN, i, j);

      return J;
    }

    FunctionReturnType evaluateEmbeddingFunctionImpl(const AnsatzFunctionType& N) const {
      FunctionReturnType res;
      setZero(res);
      for (size_t i = 0; i < coeffs.size(); ++i)
        for (size_t k = 0; k < valueSize; ++k)
          res[k] += coeffs[i].getValue()[k] * N[i];
      return res;
    }

    mutable AnsatzFunctionJacobian dNTransformed;
    Dune::CachedLocalBasis<DuneBasis> basis_;
    CoeffContainer coeffs;
    std::shared_ptr<const Geometry> geometry_;
    //    const decltype(Dune::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename DuneBasis, typename CoeffContainer, typename Geometry, std::size_t ID, typename LinAlg>
  struct LocalFunctionTraits<EmbeddedLocalFunction<DuneBasis, CoeffContainer, Geometry, ID, LinAlg>> {
    /** \brief Type used for coordinates */
    using ctype = typename CoeffContainer::value_type::ctype;
    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = CoeffContainer::value_type::valueSize;
    /** \brief Dimension of the correction size of coeffs */
    static constexpr int correctionSize = CoeffContainer::value_type::correctionSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Dune::CachedLocalBasis<DuneBasis>::gridDim;
    /** \brief The manifold where the function values lives in */
    using Manifold = typename CoeffContainer::value_type;
    /** \brief Type for the return value */
    using FunctionReturnType = typename Manifold::CoordinateType;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = typename DefaultLinearAlgebra::template FixedSizedMatrix<ctype, valueSize, gridDim>;
    /** \brief Type for the derivatives wrt. the coefficients */
    using CoeffDerivMatrix =
        typename DefaultLinearAlgebra::template FixedSizedMatrix<ctype, correctionSize, correctionSize>;
    /** \brief Type for the derivatives wrt. the coefficients */
    using CoeffDerivEukMatrix = typename DefaultLinearAlgebra::template FixedSizedMatrix<ctype, valueSize, valueSize>;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Dune::CachedLocalBasis<DuneBasis>::JacobianType;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Dune::CachedLocalBasis<DuneBasis>::AnsatzFunctionType;
    /** \brief Type for the points for evaluation, usually the integration points */
    using DomainType = typename DuneBasis::Traits::DomainType;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename DefaultLinearAlgebra::template FixedSizedVector<ctype, valueSize>;
    /** \brief Type for the directional derivatives */
    using AlongType = typename DefaultLinearAlgebra::template FixedSizedVector<ctype, valueSize>;
    /** \brief Dimension of the world where this function is mapped to from the reference element */
    static constexpr int worldDimension = Geometry::coorddimension;
  };

}  // namespace Dune
