//
// Created by Alex on 21.04.2021.
//

#pragma once

#include <dune/localfefunctions/impl/clonableLocalFunction.hh>

#include <concepts>

#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/localFunctionHelper.hh>
#include <dune/localfefunctions/localFunctionInterface.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/finiteElements/unitvectormod.hh>

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace Dune {

  template <typename DuneBasis_, typename CoeffContainer_, typename Geometry_, std::size_t ID = 0,
            typename LinAlg = Dune::DefaultLinearAlgebra>
  class GeodesicLocalFunction
      : public LocalFunctionInterface<GeodesicLocalFunction<DuneBasis_, CoeffContainer_, Geometry_, ID, LinAlg>>,
        public ClonableLocalFunction<GeodesicLocalFunction<DuneBasis_, CoeffContainer_, Geometry_, ID, LinAlg>> {
    using Interface = LocalFunctionInterface<GeodesicLocalFunction<DuneBasis_, CoeffContainer_, Geometry_, ID, LinAlg>>;

    template <size_t ID_ = 0>
    static constexpr int orderID = ID_ == ID ? 1000 : 0;

  public:
    using DuneBasis =DuneBasis_;
    using CoeffContainer =CoeffContainer_;
    using Geometry =Geometry_;
    using LinearAlgebra                                     = LinAlg;
    static constexpr bool providesDerivativeTransformations = true;
    friend Interface;
    friend ClonableLocalFunction<GeodesicLocalFunction>;
    constexpr GeodesicLocalFunction(
        const Dune::CachedLocalBasis<DuneBasis>& p_basis, const CoeffContainer& coeffs_,
        const std::shared_ptr<const Geometry>& geo,
        Dune::template index_constant<ID> = Dune::template index_constant<std::size_t(0)>{})
        : basis_{p_basis},
          coeffs{coeffs_},
          geometry_{geo} {
      coeffsReal.resize(coeffs.size());
      coeffsDual2nd.resize(coeffs.size());
      for (int i = 0; i < coeffs.size(); ++i) {
        coeffsReal[i]=coeffs[i];
        coeffsDual2nd[i]=coeffs[i];
      }
      std::cout<<"coeffsReal Constructor"<<std::endl;
      std::cout<<coeffsReal<<std::endl;
    }

    using Traits = LocalFunctionTraits<GeodesicLocalFunction>;

    static constexpr bool isLeaf = true;
    static constexpr std::array<int, 1> id{ID};




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
    struct Rebind {
      using other = GeodesicLocalFunction<
          DuneBasis, typename Std::Rebind<CoeffContainer, typename Manifold::template Rebind<OtherType>::other>::other,
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
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    static auto tryToCallThirdDerivativeOfProjectionWRTposition(const FunctionReturnType& valE, const AlongType& along,
                                                                const Eigen::Ref<const AlongType>& along2) {
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
      return  evaluateDerivativeWRTSpaceAllImplImpl(ipIndexOrPosition,transArgs,coeffs);
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    auto evaluateDerivativeWRTSpaceAllImplImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const On<TransformArgs...>& transArgs,const auto& localCoeffs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      const auto t = evaluateFunctionImplImpl(ipIndexOrPosition,transArgs,localCoeffs);
      const auto tManifold = std::remove_cvref_t<decltype(localCoeffs[0])>(t);
      Eigen::Matrix<typename std::remove_cvref_t<decltype(localCoeffs[0])>::ctype,valueSize,gridDim> gA;
      for (int i = 0; i < gridDim; ++i) {
        gA.col(i) = getGradientOfSquaredDistanceImpl(dNTransformed.col(i),localCoeffs)(tManifold);
      }
      auto H = getHessianOfSquaredDistanceImpl(N,localCoeffs)(tManifold);

      auto J              = (-H.completeOrthogonalDecomposition().pseudoInverse()*gA).eval();
      return  J;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                         int spaceIndex, const On<TransformArgs...>& transArgs) const {
      return evaluateDerivativeWRTSpaceSingleImplImpl(ipIndexOrPosition,spaceIndex,transArgs,coeffs);
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    auto evaluateDerivativeWRTSpaceSingleImplImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                         int spaceIndex, const On<TransformArgs...>& transArgs,const auto& localCoeffs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs, geometry_, ipIndexOrPosition, basis_);
      const auto t = evaluateFunctionImplImpl(ipIndexOrPosition,transArgs,localCoeffs);
      const auto tManifold = std::remove_cvref_t<decltype(localCoeffs[0])>(t);
      Eigen::Vector<typename std::remove_cvref_t<decltype(localCoeffs[0])>::ctype,valueSize> gA = getGradientOfSquaredDistanceImpl(dNTransformed.col(spaceIndex),localCoeffs)(tManifold);
      const auto H = getHessianOfSquaredDistanceImpl(N,localCoeffs)(tManifold);

      return (-H.completeOrthogonalDecomposition().pseudoInverse()*gA).eval();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           int coeffsIndex, const On<TransformArgs...>&) const {

      auto f = [&](auto& x)
      {
        coeffsReal[coeffsIndex].addInEmbedding(x);
        auto res= evaluateFunctionImplImpl(ipIndexOrPosition,On<TransformArgs...>(),coeffsReal);
        coeffsReal[coeffsIndex].addInEmbedding(-x);
        return res;
      };
      Eigen::Vector<autodiff::real,valueSize> F;
      Eigen::Vector<autodiff::real,valueSize> x;                           // the input vector x with 5 variables
      x.setZero();
      Eigen::Matrix<ctype ,valueSize,valueSize> J;
          jacobian(f, autodiff::wrt(x), autodiff::at(x), F,J);
          std::cout<<"evaluateDerivativeWRTCoeffsImpl coeffIndex:"<<coeffsIndex<<std::endl;
          std::cout<<"F"<<std::endl;
          std::cout<<F<<std::endl;
          std::cout<<"J"<<std::endl;
          std::cout<<J<<std::endl;
      return J*coeffs[coeffsIndex].orthonormalFrame();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           const std::array<size_t, 2>& coeffsIndex,
                                                           const Along<AlongArgs...>& alongArgs,
                                                           const On<TransformArgs...>& transArgs) const {
      if (coeffsIndex[0] == coeffsIndex[1])
      {
        auto f = [&](auto& x) {
          coeffsDual2nd[coeffsIndex[0]].addInEmbedding(x.template segment<valueSize>(0));
          auto res = evaluateFunctionImplImpl(ipIndexOrPosition, On<TransformArgs...>(), coeffsDual2nd);
          coeffsDual2nd[coeffsIndex[0]].addInEmbedding(-x.template segment<valueSize>(0));
          return res.dot(std::get<0>(alongArgs.args));
        };
        Eigen::Vector<ctype, valueSize> g;
        Eigen::Matrix<ctype, valueSize,valueSize> H;
        autodiff::dual2nd e;
        Eigen::Vector<autodiff::dual2nd, valueSize> x;  // the input vector x with 5 variables
        x.setZero();
        hessian(f, autodiff::wrt(x), autodiff::at(x), e, g, H);
        if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
          H -= (coeffs[coeffsIndex[0]].getValue().dot(g)) * CoeffDerivEukMatrix::Identity();
        }
        std::cout<<"evaluateSecondDerivativeWRTCoeffsImpl coeffIndex:"<<coeffsIndex[0]<<" "<<coeffsIndex[1]<<std::endl;
        std::cout<<"g"<<std::endl;
        std::cout<<g<<std::endl;
        std::cout<<"H"<<std::endl;
        std::cout<<H<<std::endl;
        return coeffs[coeffsIndex[0]].orthonormalFrame().transpose() * H * coeffs[coeffsIndex[0]].orthonormalFrame();
      }else

      {
        auto f = [&](auto& x) {
          coeffsDual2nd[coeffsIndex[0]].addInEmbedding(x.template segment<valueSize>(0));
            coeffsDual2nd[coeffsIndex[1]].addInEmbedding(x.template segment<valueSize>(valueSize));
          auto res = evaluateFunctionImplImpl(ipIndexOrPosition, On<TransformArgs...>(), coeffsDual2nd);
          coeffsDual2nd[coeffsIndex[0]].addInEmbedding(-x.template segment<valueSize>(0));
            coeffsDual2nd[coeffsIndex[1]].addInEmbedding(-x.template segment<valueSize>(valueSize));
          return res.dot(std::get<0>(alongArgs.args));
        };
        Eigen::Vector<ctype, 2 * valueSize> g;
        Eigen::Matrix<ctype, 2 * valueSize, 2 * valueSize> H;
        autodiff::dual2nd e;
        Eigen::Vector<autodiff::dual2nd, 2 * valueSize> x;  // the input vector x with 5 variables
        x.setZero();
        hessian(f, autodiff::wrt(x), autodiff::at(x), e, g, H);
        CoeffDerivEukMatrix Hr = H.template block<valueSize, valueSize>(0, valueSize);
        std::cout<<"evaluateSecondDerivativeWRTCoeffsImpl diff coeffIndex:"<<coeffsIndex[0]<<" "<<coeffsIndex[1]<<std::endl;
        std::cout<<"g"<<std::endl;
        std::cout<<g<<std::endl;
        std::cout<<"H"<<std::endl;
        std::cout<<H<<std::endl;
        return coeffs[coeffsIndex[0]].orthonormalFrame().transpose() * Hr * coeffs[coeffsIndex[1]].orthonormalFrame();
      }
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    std::array<CoeffDerivEukRieMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const On<TransformArgs...>& transArgs) const {
      //      std::array<CoeffDerivEukMatrix, gridDim> WarrayEuk
      //          = evaluateDerivativeWRTCoeffsANDSpatialEukImpl(ipIndexOrPosition, coeffsIndex, transArgs);
      std::array<CoeffDerivEukRieMatrix, gridDim> WarrayRie;
      //      const auto BLA = coeffs[coeffsIndex].orthonormalFrame().eval();
      //      for (int dir = 0; dir < gridDim; ++dir)
      //        WarrayRie[dir] = WarrayEuk[dir] * BLA;
      std::cout<<"coeffsReal Outside"<<std::endl;
      std::cout<<coeffsReal<<std::endl;
      auto f = [&](auto& x)
      {
        std::cout<<"coeffsReal Inside Before"<<std::endl;
        std::cout<<coeffsReal<<std::endl;
        std::cout<<"x"<<std::endl;
        std::cout<<x<<std::endl;
        coeffsReal[coeffsIndex].addInEmbedding(x);
        auto J = evaluateDerivativeWRTSpaceAllImplImpl(ipIndexOrPosition,On<TransformArgs...>(),coeffsReal);
        coeffsReal[coeffsIndex].addInEmbedding(-x);
        std::cout<<"coeffsReal Inside After"<<std::endl;
        std::cout<<coeffsReal<<std::endl;
        auto JV = J.reshaped().eval();
        return JV;
      };
      Eigen::Vector<autodiff::real,2*valueSize> F;
      Eigen::Vector<autodiff::real,valueSize> x;
      x.setZero();
      Eigen::Matrix<ctype ,2*valueSize,valueSize> JV;
      jacobian(f, autodiff::wrt(x), autodiff::at(x), F,JV);
      for (int i = 0; i < gridDim; ++i) {
        WarrayRie[i]=JV.template block<valueSize,valueSize>(valueSize*i,0)*coeffs[coeffsIndex].orthonormalFrame();
      }
      std::cout<<"evaluateDerivativeWRTCoeffsANDSpatialImpl coeffIndex:"<<coeffsIndex<<std::endl;
      std::cout<<"JV"<<std::endl;
      std::cout<<JV<<std::endl;
      std::cout<<"WarrayRie[0]"<<std::endl;
      std::cout<<WarrayRie[0]<<std::endl;
      std::cout<<"WarrayRie[1]"<<std::endl;
      std::cout<<WarrayRie[1]<<std::endl;
      return WarrayRie;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
        const On<TransformArgs...>& transArgs) const
      {
        CoeffDerivEukRieMatrix WarrayRie;

        auto f = [&](auto& x) {
          coeffsReal[coeffsIndex].addInEmbedding(x);
          auto J  = evaluateDerivativeWRTSpaceSingleImplImpl(ipIndexOrPosition,spatialIndex, On<TransformArgs...>(), coeffsReal);
          coeffsReal[coeffsIndex].addInEmbedding(-x);

          return J;
        };
        Eigen::Vector<autodiff::real, valueSize> F;
        Eigen::Vector<autodiff::real, valueSize> x;  // the input vector x with 5 variables
        x.setZero();
        Eigen::Matrix<ctype, valueSize, valueSize> JV;
        jacobian(f, wrt(x), at(x), F, JV);
          WarrayRie
              = JV* coeffs[coeffsIndex].orthonormalFrame();
          std::cout<<"evaluateDerivativeWRTCoeffsANDSpatialSingleImpl coeffIndex:"<<coeffsIndex<<"spatialindex"<<spatialIndex<<std::endl;
          std::cout<<"JV"<<std::endl;
          std::cout<<JV<<std::endl;
          std::cout<<"WarrayRi"<<std::endl;
          std::cout<<WarrayRie<<std::endl;
        return WarrayRie;
    }

//    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
//    CoeffDerivEukMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(
//        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
//        const On<TransformArgs...>& transArgs) const {
//      //      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
//      //      maytransformDerivatives(dNraw, dNTransformed, transArgs);
//      //      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
//      //      const JacobianColType Jcol    = evaluateEmbeddingJacobianColImpl(dNTransformed, spatialIndex);
//      //      const CoeffDerivEukMatrix Pm  = tryToCallDerivativeOfProjectionWRTposition(valE);
//      CoeffDerivEukMatrix W;
//      //      const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, Jcol);
//      //      W             = Qi * N[coeffsIndex] + Pm * dNTransformed(coeffsIndex, spatialIndex);
//      return W;
//    }
    template <bool isSingle,typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialGenericImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const Along<AlongArgs...>& alongArgs, const On<TransformArgs...>& transArgs, std::optional<int> spaceIndex) const {
        const auto& along             = std::get<0>(alongArgs.args);
        if (coeffsIndex[0] == coeffsIndex[1]) {
        auto f = [&](auto& x)
        {
          coeffsDual2nd[coeffsIndex[0]].addInEmbedding(x.template segment<valueSize>(0));
          if constexpr (isSingle)
          {
            assert(spaceIndex.has_value());
            auto res = evaluateDerivativeWRTSpaceSingleImplImpl(ipIndexOrPosition,spaceIndex.value(), On<TransformArgs...>(), coeffsDual2nd);
            coeffsDual2nd[coeffsIndex[0]].addInEmbedding(-x.template segment<valueSize>(0));
            return (res.transpose() * along).trace();
          }else {
            auto res = evaluateDerivativeWRTSpaceAllImplImpl(ipIndexOrPosition, On<TransformArgs...>(), coeffsDual2nd);
            coeffsDual2nd[coeffsIndex[0]].addInEmbedding(-x.template segment<valueSize>(0));
            return (res.transpose() * along).trace();
          }
        };
        Eigen::Vector<ctype, valueSize> g;
        Eigen::Matrix<ctype, valueSize, valueSize> H;
        autodiff::dual2nd e;
        Eigen::Vector<autodiff::dual2nd, valueSize> x;  // the input vector x with 5 variables
        x.setZero();
        hessian(f, autodiff::wrt(x), autodiff::at(x), e, g, H);


        // Riemannian Hessian Weingarten map correction
        H      -= coeffs[coeffsIndex[0]].getValue().dot(g) * CoeffDerivEukMatrix::Identity();
        std::cout<<"evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialGenericImpl coeffIndex:"<<coeffsIndex[0]<<" "<<coeffsIndex[1]<<"spatialindex"<<spaceIndex.value()<<std::endl;
        std::cout<<"H"<<std::endl;
        std::cout<<H<<std::endl;
        std::cout<<"g"<<std::endl;
        std::cout<<g<<std::endl;
        CoeffDerivMatrix ChiArrayRie;
        const auto BLA0 = coeffs[coeffsIndex[0]].orthonormalFrame();
        ChiArrayRie      = BLA0.transpose() * H * BLA0;
        return ChiArrayRie;
        }else
        {

        auto f = [&](auto& x)
        {
          coeffsDual2nd[coeffsIndex[0]].addInEmbedding(x.template segment<valueSize>(0));
          coeffsDual2nd[coeffsIndex[1]].addInEmbedding(x.template segment<valueSize>(valueSize));
          auto res= evaluateDerivativeWRTSpaceAllImplImpl(ipIndexOrPosition,On<TransformArgs...>(),coeffsDual2nd);
          coeffsDual2nd[coeffsIndex[0]].addInEmbedding(-x.template segment<valueSize>(0));
          coeffsDual2nd[coeffsIndex[1]].addInEmbedding(-x.template segment<valueSize>(valueSize));
          return (res.transpose()*along).trace();
        };
        Eigen::Vector<ctype, 2 * valueSize> g;
        Eigen::Matrix<ctype, 2 * valueSize, 2 * valueSize> H;
        autodiff::dual2nd e;
        Eigen::Vector<autodiff::dual2nd, 2 * valueSize> x;  // the input vector x with 5 variables
        x.setZero();
        hessian(f, autodiff::wrt(x), autodiff::at(x), e, g, H);
        CoeffDerivEukMatrix Hr = H.template block<valueSize, valueSize>(0, valueSize);
        std::cout<<"evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialGenericImpl coeffIndex:"<<coeffsIndex[0]<<" "<<coeffsIndex[1]<<std::endl;
        std::cout<<"H"<<std::endl;
        std::cout<<H<<std::endl;
        std::cout<<"Hr"<<std::endl;
        std::cout<<Hr<<std::endl;
        std::cout<<"g"<<std::endl;
        std::cout<<g<<std::endl;
        CoeffDerivMatrix ChiArrayRie;
        const auto BLA0T = coeffs[coeffsIndex[0]].orthonormalFrame().transpose().eval();
        const auto BLA1  = coeffs[coeffsIndex[1]].orthonormalFrame();
        ChiArrayRie      = BLA0T * Hr * BLA1;
        return ChiArrayRie;
        }
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const Along<AlongArgs...>& alongArgs, const On<TransformArgs...>& transArgs) const {
        return evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialGenericImpl<false>(ipIndexOrPosition,coeffsIndex,alongArgs,transArgs);
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const int spatialIndex, const Along<AlongArgs...>& alongArgs, const On<TransformArgs...>& transArgs) const {
        return evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialGenericImpl<true>(ipIndexOrPosition,coeffsIndex,alongArgs,transArgs,spatialIndex);
    }


    auto getGradientOfSquaredDistance(const auto& N) const
    {
      return getGradientOfSquaredDistanceImpl(N,coeffs);
    }

    auto getGradientOfSquaredDistanceImpl(const auto& N, const auto& localCoeffs) const
    {
      auto grad= [&](auto& t){
        Eigen::Vector<typename std::remove_cvref_t<decltype(localCoeffs[0])>::ctype,valueSize> gradVec;
        gradVec.setZero();
        for (size_t i = 0; i < localCoeffs.size(); ++i)
          gradVec+=euclideanDerivativeOfDistanceSquaredWRTSecondArgument(localCoeffs[i],t)*N[i];

        return gradVec;
      };
      return grad;
    }

    auto getHessianOfSquaredDistance(const auto& N) const
    {
      return getHessianOfSquaredDistanceImpl(N,coeffs);
    }

    auto getHessianOfSquaredDistanceImpl(const auto& N, const auto& localCoeffs) const
    {
      auto hess= [&](auto& t){
        Eigen::Matrix<typename std::remove_cvref_t<decltype(localCoeffs[0])>::ctype,valueSize,valueSize> hessM;
        hessM.setZero();
        for (size_t i = 0; i < localCoeffs.size(); ++i) {
          hessM+=euclideanSecondDerivativeOfDistanceSquaredWRTSecondArgument(localCoeffs[i],t)*N[i];
        }
        return hessM;
      };
      return hess;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    FunctionReturnType evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                            [[maybe_unused]] const On<TransformArgs...>&) const {
      return evaluateFunctionImplImpl(ipIndexOrPosition,On<TransformArgs...>(),coeffs);
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    auto evaluateFunctionImplImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                            [[maybe_unused]] const On<TransformArgs...>&, const auto& localCoeffs) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      using ManifoldLocal= std::remove_cvref_t<decltype(localCoeffs[0])>;
      typename ManifoldLocal::CoordinateType tresV;
      tresV.setZero();
//      =localCoeffs[0];

      for (int i = 0; i < localCoeffs.size(); ++i) {
        tresV+=N[i]*localCoeffs[i].getValue();
      }
      auto tres= ManifoldLocal(tresV);
      // predictor
      auto g = getGradientOfSquaredDistanceImpl(N,localCoeffs);
      auto h = getHessianOfSquaredDistanceImpl(N,localCoeffs);
      int maxiter=20;
      int miniter=0;
      int i=0;
      typename decltype(tres)::ctype gRNorm;
      for (; i<maxiter; ++i)
      {
        auto gE=g(tres);
        auto P = tres.derivativeOfProjectionWRTposition(tres.getValue());
        auto hE=(P*h(tres)*P).eval();
        auto d=(-hE.completeOrthogonalDecomposition().pseudoInverse() * gE).eval();
        tres.setValue(tres.getValue()+P*d);

        gRNorm=(P*gE).norm();
        if(gRNorm<1e-8 and i>= miniter)
          break;
      }
//      auto nonLinOp  = Ikarus::NonLinearOperator( Ikarus::functions(g,h), Ikarus::parameter(tres));
//      Ikarus::NewtonRaphson nr(nonLinOp, [&](auto& r, auto& A_) { return A_.completeOrthogonalDecomposition().pseudoInverse() * r; },
//          [&](auto& x, const auto& b) { x.setValue(x.getValue()+b); });
//      nr.setup({.Rtol=1e-8});
//      const auto solverInfo =  nr.solve();
      if(i==maxiter-1)
        DUNE_THROW(Dune::InvalidStateException,"GFE interpolation not converged! Residual norm is "<<gRNorm);
      return tres.getValue();
    }


    mutable AnsatzFunctionJacobian dNTransformed;
    Dune::CachedLocalBasis<DuneBasis> basis_;
    CoeffContainer coeffs;
    mutable Dune::BlockVector<typename CoeffContainer::value_type::template rebind<autodiff::real >::other> coeffsReal;
    mutable Dune::BlockVector<typename CoeffContainer::value_type::template rebind<autodiff::dual2nd >::other> coeffsDual2nd;
    std::shared_ptr<const Geometry> geometry_;
  };

  template <typename DuneBasis, typename CoeffContainer, typename Geometry, std::size_t ID,
            typename LinAlg>
  struct LocalFunctionTraits<GeodesicLocalFunction<DuneBasis, CoeffContainer,Geometry, ID,LinAlg>> {
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

}  // namespace Ikarus
