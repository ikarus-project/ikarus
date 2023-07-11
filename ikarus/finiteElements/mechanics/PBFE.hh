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
  class ProjectionBasedLocalFunction2
      : public LocalFunctionInterface<ProjectionBasedLocalFunction2<DuneBasis_, CoeffContainer_, Geometry_, ID, LinAlg>>,
        public ClonableLocalFunction<ProjectionBasedLocalFunction2<DuneBasis_, CoeffContainer_, Geometry_, ID, LinAlg>> {
    using Interface = LocalFunctionInterface<ProjectionBasedLocalFunction2>;

    template <size_t ID_ = 0>
    static constexpr int orderID = ID_ == ID ? 1000 : 0;

  public:
    using DuneBasis =DuneBasis_;
    using CoeffContainer =CoeffContainer_;
    using Geometry =Geometry_;
    friend Interface;
    friend ClonableLocalFunction<ProjectionBasedLocalFunction2>;
    constexpr ProjectionBasedLocalFunction2(
        const Dune::CachedLocalBasis<DuneBasis>& p_basis, const CoeffContainer& coeffs_,
        const std::shared_ptr<const Geometry>& geo,
        Dune::template index_constant<ID> = Dune::template index_constant<std::size_t(0)>{})
        : basis_{p_basis},
          coeffs{coeffs_},
          geometry_{geo}  //          ,coeffsAsMat{Dune::viewAsEigenMatrixFixedDyn(coeffs)}
    {}

    using Traits = LocalFunctionTraits<ProjectionBasedLocalFunction2>;

    static constexpr bool isLeaf = true;
    static constexpr std::array<int, 1> id{ID};

    using LinearAlgebra = LinAlg;

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
      using other = ProjectionBasedLocalFunction2<
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

    template <typename DomainTypeOrIntegrationPointIndex, typename TransformArgs>
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const On<TransformArgs>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      return tryToCallDerivativeOfProjectionWRTposition(valE) * J;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename TransformArgs>
    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                         int spaceIndex, const On<TransformArgs>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      JacobianColType Jcol    = evaluateEmbeddingJacobianColImpl(dNTransformed, spaceIndex);
      FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      return tryToCallDerivativeOfProjectionWRTposition(valE) * Jcol;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename TransformArgs>
    CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           int coeffsIndex, const On<TransformArgs>& transArgs) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      return (evaluateDerivativeWRTCoeffsEukImpl(N, coeffsIndex) * coeffs[coeffsIndex].orthonormalFrame());
    }

    CoeffDerivEukMatrix evaluateDerivativeWRTCoeffsEukImpl(const AnsatzFunctionType& N, int coeffsIndex) const {
      FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      return tryToCallDerivativeOfProjectionWRTposition(valE) * N[coeffsIndex];
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename TransformArgs>
    CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           const std::array<size_t, 2>& coeffsIndex,
                                                           const Along<AlongArgs...>& alongArgs,
                                                           const On<TransformArgs>& transArgs) const {
      const auto& N                 = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);

      CoeffDerivEukMatrix ddt = tryToCallSecondDerivativeOfProjectionWRTposition(valE, std::get<0>(alongArgs.args))
                                * N[coeffsIndex[0]] * N[coeffsIndex[1]];

      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const CoeffDerivEukMatrix dt = evaluateDerivativeWRTCoeffsEukImpl(N, coeffsIndex[0]);
        auto& unitVec                = coeffs[coeffsIndex[0]];
        auto& unitVecVal             = unitVec.getValue();
        auto scal                    = inner(unitVecVal, dt * std::get<0>(alongArgs.args));
        auto idmat                   = createScaledIdentityMatrix<ctype, valueSize, valueSize>(scal);
        ddt -= idmat;
      }

      return transposeEvaluated(coeffs[coeffsIndex[0]].orthonormalFrame()) * ddt
             * coeffs[coeffsIndex[1]].orthonormalFrame();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename TransformArgs>
    std::array<CoeffDerivEukRieMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const On<TransformArgs>& transArgs) const {
      std::array<CoeffDerivEukMatrix, gridDim> WarrayEuk
          = evaluateDerivativeWRTCoeffsANDSpatialEukImpl(ipIndexOrPosition, coeffsIndex, transArgs);
      std::array<CoeffDerivEukRieMatrix, gridDim> WarrayRie;
      const auto BLA = coeffs[coeffsIndex].orthonormalFrame();
      for (int dir = 0; dir < gridDim; ++dir)
        WarrayRie[dir] = WarrayEuk[dir] * BLA;

      return WarrayRie;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename TransformArgs>
    std::array<CoeffDerivEukMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialEukImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const On<TransformArgs>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      const CoeffDerivEukMatrix Pm  = tryToCallDerivativeOfProjectionWRTposition(valE);
      std::array<CoeffDerivEukMatrix, gridDim> Warray;
      for (int dir = 0; dir < gridDim; ++dir) {
        const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, col(J, dir));
        Warray[dir]   = Qi * N[coeffsIndex] + Pm * coeff(dNTransformed, coeffsIndex, dir);
      }

      return Warray;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename TransformArgs>
    CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
        const On<TransformArgs>& transArgs) const {
      const CoeffDerivEukMatrix WEuk
          = evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(ipIndexOrPosition, coeffsIndex, spatialIndex, transArgs);
      return WEuk * coeffs[coeffsIndex].orthonormalFrame();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename TransformArgs>
    CoeffDerivEukMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
        const On<TransformArgs>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      const JacobianColType Jcol    = evaluateEmbeddingJacobianColImpl(dNTransformed, spatialIndex);
      const CoeffDerivEukMatrix Pm  = tryToCallDerivativeOfProjectionWRTposition(valE);
      CoeffDerivEukMatrix W;
      const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, Jcol);
      W             = Qi * N[coeffsIndex] + Pm * coeff(dNTransformed, coeffsIndex, spatialIndex);
      return W;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename TransformArgs>
    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const Along<AlongArgs...>& alongArgs, const On<TransformArgs>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);

      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      const auto& along             = std::get<0>(alongArgs.args);
      CoeffDerivEukMatrix ChiArrayEuk;
      setZero(ChiArrayEuk);
      static_assert(Cols<decltype(along)>::value == gridDim);
      static_assert(Cols<decltype(J)>::value == gridDim);
      static_assert(Rows<decltype(along)>::value == valueSize);
      static_assert(Rows<decltype(J)>::value == valueSize);
      for (int i = 0; i < gridDim; ++i) {
        auto colAlong               = col(along, i);
        const auto chi              = tryToCallThirdDerivativeOfProjectionWRTposition(valE, colAlong, col(J, i));
        const CoeffDerivEukMatrix S = tryToCallSecondDerivativeOfProjectionWRTposition(valE, colAlong);
        const auto& NI              = N[coeffsIndex[0]];
        const auto& NJ              = N[coeffsIndex[1]];
        const auto& dNIdi           = coeff(dNTransformed, coeffsIndex[0], i);
        const auto& dNJdi           = coeff(dNTransformed, coeffsIndex[1], i);
        ChiArrayEuk += chi * NI * NJ;
        ChiArrayEuk += S * (dNIdi * NJ + dNJdi * NI);
      }
      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const std::array<CoeffDerivEukMatrix, gridDim> Warray
            = evaluateDerivativeWRTCoeffsANDSpatialEukImpl(ipIndexOrPosition, coeffsIndex[0], transArgs);
        for (int i = 0; i < gridDim; ++i) {
          ChiArrayEuk -= createScaledIdentityMatrix<ctype, valueSize, valueSize>(
              inner(coeffs[coeffsIndex[0]].getValue(), Warray[i] * col(along, i)));
        }
      }

      CoeffDerivMatrix ChiArrayRie;
      const auto BLA0T = eval(transposeEvaluated(coeffs[coeffsIndex[0]].orthonormalFrame()));
      const auto BLA1  = coeffs[coeffsIndex[1]].orthonormalFrame();
      ChiArrayRie      = BLA0T * ChiArrayEuk * BLA1;
      return ChiArrayRie;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename TransformArgs>
    CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const int spatialIndex, const Along<AlongArgs...>& alongArgs, const On<TransformArgs>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      const auto& along             = std::get<0>(alongArgs.args);
      const CoeffDerivEukMatrix S   = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along);
      const auto chi          = tryToCallThirdDerivativeOfProjectionWRTposition(valE, along, col(J, spatialIndex));
      const auto& NI          = N[coeffsIndex[0]];
      const auto& NJ          = N[coeffsIndex[1]];
      const auto& dNIdi       = coeff(dNTransformed, coeffsIndex[0], spatialIndex);
      const auto& dNJdi       = coeff(dNTransformed, coeffsIndex[1], spatialIndex);
      CoeffDerivEukMatrix Chi = chi * NI * NJ;
      Chi += S * (dNIdi * NJ + dNJdi * NI);

      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const CoeffDerivEukMatrix W = evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(
            ipIndexOrPosition, coeffsIndex[0], spatialIndex, transArgs);
        Chi -= createScaledIdentityMatrix<ctype, valueSize, valueSize>(
            inner(coeffs[coeffsIndex[0]].getValue(), W * along));
      }
      return (transposeEvaluated(coeffs[coeffsIndex[0]].orthonormalFrame()) * Chi
              * coeffs[coeffsIndex[1]].orthonormalFrame());
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename TransformArgs>
    FunctionReturnType evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                            [[maybe_unused]] const On<TransformArgs>&) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      return Manifold(evaluateEmbeddingFunctionImpl(N)).getValue();
    }

    JacobianColType evaluateEmbeddingJacobianColImpl(const AnsatzFunctionJacobian& dN, int spaceIndex) const {
      JacobianColType Jcol;
      setZero(Jcol);
      for (int j = 0; j < Rows<JacobianColType>::value; ++j) {
        for (int i = 0; i < coeffs.size(); ++i)
          Jcol[j] += coeffs[i].getValue()[j] * coeff(dN, i, spaceIndex);
      }

      return Jcol;
    }

    Jacobian evaluateEmbeddingJacobianImpl(const AnsatzFunctionJacobian& dN) const {
      Jacobian J;
      setZero(J);
      for (int j = 0; j < gridDim; ++j)
        for (int k = 0; k < valueSize; ++k)
          for (int i = 0; i < coeffs.size(); ++i)
            coeff(J, k, j) += coeffs[i].getValue()[k] * coeff(dN, i, j);

      return J;
    }

    FunctionReturnType evaluateEmbeddingFunctionImpl(const AnsatzFunctionType& N) const {
      FunctionReturnType res;
      setZero(res);
      for (int i = 0; i < coeffs.size(); ++i)
        for (int k = 0; k < valueSize; ++k)
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
  struct LocalFunctionTraits<ProjectionBasedLocalFunction2<DuneBasis, CoeffContainer, Geometry, ID, LinAlg>> {
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
