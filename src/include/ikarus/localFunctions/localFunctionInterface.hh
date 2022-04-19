//
// Created by alex on 3/17/22.
//

#pragma once
#include <concepts>

#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/utils/traits.hh>
namespace Ikarus {
  namespace DerivativeDirections {
  static struct Spatialall{
    using value = std::integral_constant<int, -2>;
  } spatialall;

    struct SpatialPartial {
      using value = std::integral_constant<int, -3>;
      size_t index{};
    };

  struct SingleCoeff {
    using value = std::integral_constant<int, -4>;
    size_t index{};
  };

  struct TwoCoeff {
    using value = std::integral_constant<int, -5>;
    std::array<size_t,2> index{};
  };

    SpatialPartial spatial(size_t i);
     SingleCoeff coeff(size_t i);
  TwoCoeff coeff(size_t i,size_t j);

    template <int Dim>
    struct Counter {
      int coeffDerivatives{};
      int spatialall{};
      int dynamicspatial{};
    };

  template <int gridDim, int coeffDerivs=2>
  struct NewCounter {
    std::array<int, gridDim> spatialall{};
    int dynamicspatial{};
  };


    template <typename WrtType, int gridDim>
    std::array<int, gridDim> countDynamicSpatialDerivativesInTuple(const WrtType& wrt) {
      std::array<int, gridDim> counter{};
      Dune::Hybrid::forEach(
          Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(WrtType::args)>>()), [&](auto i) {
            using currentDerivType = std::tuple_element_t<i, typename WrtType::Args>;
            if constexpr (std::is_same_v<currentDerivType, SpatialPartial>) ++counter[std::get<i>(wrt.args).index];
          });
      return counter;
    }


    struct ConstExprCounter
    {
        int singleCoeffDerivs{};
        int twoCoeffDerivs{};
        int spatialDerivs{};
        int spatialAll{};
    };

  template <typename WrtType, int gridDim>
  consteval ConstExprCounter countDerivativesType() {

    ConstExprCounter counter{};
    using Tuple = typename WrtType::Args;
    counter.singleCoeffDerivs = Ikarus::Std::countType<Tuple,SingleCoeff>();
    counter.twoCoeffDerivs = Ikarus::Std::countType<Tuple,TwoCoeff>();
    counter.spatialDerivs = Ikarus::Std::countType<Tuple,SpatialPartial>();
    counter.spatialAll = Ikarus::Std::countType<Tuple,Spatialall>();
    return counter;
  }

  template <typename WrtType>
  auto extractCoeffIndices(WrtType&& wrt) {

    if constexpr(Ikarus::Std::hasType<SingleCoeff,typename WrtType::Args>::value)
    return std::get<SingleCoeff>(wrt.args).index; //returns single int
    else if constexpr(Ikarus::Std::hasType<TwoCoeff,typename WrtType::Args>::value)
    return std::get<TwoCoeff>(wrt.args).index; //return std::array<size_t,2>
  }


    template <size_t gridDim>
    int findSingleSpatial(const std::array<int, gridDim>& spatialDerivs) {
      const auto element = std::ranges::find_if_not(spatialDerivs, [](auto&& el) { return el == 0; });
      return std::distance(spatialDerivs.begin(), element);
    }

  template <typename WrtType, int gridDim>
  concept HasTwoCoeff = (countDerivativesType<WrtType, gridDim>().twoCoeffDerivs == 1);

    template <typename WrtType, int gridDim>
  concept HasSingleCoeff = (countDerivativesType<WrtType, gridDim>().singleCoeffDerivs == 1);

  template <typename WrtType, int gridDim>
  concept HasNoCoeff = (countDerivativesType<WrtType, gridDim>().singleCoeffDerivs == 0 and countDerivativesType<WrtType, gridDim>().twoCoeffDerivs == 0);

    template <typename WrtType, int gridDim>
    concept HasNoSpatial
        = (countDerivativesType<WrtType, gridDim>().spatialDerivs == 0 and countDerivativesType<WrtType, gridDim>().spatialAll == 0);

    template <typename WrtType, int gridDim>
    concept HasOneSpatialAll = countDerivativesType<WrtType, gridDim>().spatialAll == 1;

    template <typename WrtType, int gridDim>
    concept HasOneSpatialSingle = countDerivativesType<WrtType, gridDim>().spatialDerivs == 1;

    template <typename WrtType, int gridDim>
    concept HasOneSpatial = HasOneSpatialSingle<WrtType, gridDim> or HasOneSpatialAll<WrtType, gridDim>;

  }  // namespace DerivativeDirections

  template <typename... Args_>
  struct Wrt {
    using Args = std::tuple<std::remove_cvref_t<Args_>...>;
    Args args;
  };

  template <typename... Args>
  auto wrt(Args&&... args) {
    return Wrt<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
  }

  template <typename... Args>
  struct Along {
    std::tuple<Args...> args;
  };

  template <typename... Args>
  auto along(Args&&... args) {
    return Along<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
  }

  template <typename... Args>
  struct TransformWith {
    std::tuple<Args...> args;
  };

  template <typename... Args>
  auto transformWith(Args&&... args) {
    return TransformWith<Args&&...>({std::forward<Args>(args)...});
  }

  template <typename Derived>
  struct LocalFunctionTraits;

  template <typename LocalFunctionImpl>
  concept HasevaluateSecondDerivativeWRTCoeffs = requires(LocalFunctionImpl func) {
    func.evaluateSecondDerivativeWRTCoeffs(
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>(),
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AlongType>(),
        std::declval<std::array<size_t, 2>>());
  };

  template <typename LocalFunctionImpl>
  concept HasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl = requires(LocalFunctionImpl func) {
    func.evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>(),
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AlongType>(),
        std::declval<std::array<size_t, 2>>());
  };

  template <typename LocalFunctionImpl>
  concept HasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl = requires(LocalFunctionImpl func) {
    func.evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>(),
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AlongType>(),
        std::declval<std::array<size_t, 2>>(), std::declval<int>());
  };

  template <typename LocalFunctionImpl>
  concept HasevaluateDerivativeWRTSpaceAllImpl = requires(LocalFunctionImpl func) {
    func.evaluateDerivativeWRTSpaceAllImpl(
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>());
  };

  template <typename LocalFunctionImpl>
  concept HasevaluateDerivativeWRTSpaceSingleImpl = requires(LocalFunctionImpl func) {
    func.evaluateDerivativeWRTSpaceSingleImpl(
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
        std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>(), std::declval<int>());
  };

  template <typename LocalFunctionImpl>
  class LocalFunctionInterface {
  public:
    using Traits                 = LocalFunctionTraits<LocalFunctionImpl>;
    using DomainType             = typename Traits::DomainType;
    using ctype                  = typename Traits::ctype;
    using FunctionReturnType     = typename Traits::FunctionReturnType;
    using AnsatzFunctionType     = typename Traits::AnsatzFunctionType;
    using Jacobian               = typename Traits::Jacobian;
    using JacobianColType        = typename Traits::JacobianColType;
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;
    using CoeffDerivMatrix       = typename Traits::CoeffDerivMatrix;

    // Check the capabilities of derived implementation
    static constexpr bool hasSecondDerivativeWRTCoeffs = HasevaluateSecondDerivativeWRTCoeffs<LocalFunctionImpl>;
    static constexpr bool hasThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl
        = HasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl<LocalFunctionImpl>;
    static constexpr bool hasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl
        = HasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl<LocalFunctionImpl>;
    static constexpr bool hasevaluateDerivativeWRTSpaceAllImpl
        = HasevaluateDerivativeWRTSpaceAllImpl<LocalFunctionImpl>;
    static constexpr bool hasevaluateDerivativeWRTSpaceSingleImpl
        = HasevaluateDerivativeWRTSpaceSingleImpl<LocalFunctionImpl>;

    static constexpr int gridDim = Traits::gridDim;
    using TransformMatrix        = Eigen::Matrix<double, gridDim, gridDim>;

    template <typename WrtType>
    static constexpr bool hasTwoCoeff = DerivativeDirections::HasTwoCoeff<WrtType, gridDim>;
    template <typename WrtType>
    static constexpr bool hasSingleCoeff = DerivativeDirections::HasSingleCoeff<WrtType, gridDim>;
    template <typename WrtType>
    static constexpr bool hasNoCoeff = DerivativeDirections::HasNoCoeff<WrtType, gridDim>;

    /** \brief Forward the binding to the local basis */
    template <typename IntegrationRule, typename... Ints>
    requires std::conjunction_v<std::is_convertible<int, Ints>...>
    void bind(IntegrationRule&& p_rule, Impl::Derivatives<Ints...>&& ints) {
      impl().basis.bind(std::forward<IntegrationRule>(p_rule), std::forward<Impl::Derivatives<Ints...>>(ints));
    }

    /** \brief Return the function value at the i-th bound integration point*/
    FunctionReturnType evaluateFunction(long unsigned i) {
      const auto& N = impl().basis.evaluateFunction(i);
      return impl().evaluateFunctionImpl(N);
    }

    /** \brief Return the function value at the coordinates local */
    FunctionReturnType evaluateFunction(const DomainType& local) {
      AnsatzFunctionType N;
      impl().basis.evaluateFunction(local, N);
      return impl().evaluateFunctionImpl(N);
    }

    /** \brief Function to forward the call of no spatial derivative and two derivative wrt. coefficients.
     * You have to pass a along argument which specifies the direction wher this derivative is applied
     */
    template <typename... Args, typename... TransformArgs, typename... AlongArgs,
              typename DomainTypeOrIntegrationPointIndex>
    requires(
        hasTwoCoeff<Wrt<Args...>> and DerivativeDirections::HasNoSpatial<Wrt<Args...>, gridDim>and
            hasSecondDerivativeWRTCoeffs) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                                  Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                                  TransformWith<TransformArgs...>&& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);

      const auto coeffsIndices = extractCoeffIndices(std::forward<Wrt<Args...>>(args));

      // Check if a matrix is given to transform derivatives. Otherwise we do nothing
      if constexpr (sizeof...(TransformArgs) > 0) {
        AnsatzFunctionJacobian dN = (dNraw * std::get<0>(transArgs.args)).eval();
        return impl().evaluateSecondDerivativeWRTCoeffs(N, dN, std::get<0>(along.args), coeffsIndices);
      } else
        return impl().evaluateSecondDerivativeWRTCoeffs(N, dNraw, std::get<0>(along.args), coeffsIndices);
    }

    /** \brief Function to forward the call of no spatial derivative and two derivative wrt. coefficients.
     * You have to pass a along argument which specifies the direction where this derivative is applied
     * Specialization when no transformWith is passed
     */
    template <typename... Args, typename... AlongArgs, typename DomainTypeOrIntegrationPointIndex>
    requires(hasTwoCoeff<Wrt<Args...>> and DerivativeDirections::HasNoSpatial<Wrt<Args...>, gridDim>and
                 HasevaluateSecondDerivativeWRTCoeffs<
                     LocalFunctionImpl>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                                 Wrt<Args...>&& args, Along<AlongArgs...>&& along) const {

      return evaluateDerivative(localOrIpId, std::forward<Wrt<Args...>>(args), std::forward<Along<AlongArgs...>>(along),
                                transformWith());
    }

    /** \brief Function to forward the call of one spatial derivative in all directions and two derivative wrt.
     * coefficients. You have to pass a along argument which specifies the direction where this derivative is applied */
    template <typename... Args, typename... TransformArgs, typename... AlongArgs,
              typename DomainTypeOrIntegrationPointIndex>
    requires(
        hasTwoCoeff<Wrt<Args...>>and DerivativeDirections::HasOneSpatialAll<Wrt<Args...>, gridDim>and
            hasThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex&
                                                                                           localOrIpId,
                                                                                       Wrt<Args...>&& args,
                                                                                       Along<AlongArgs...>&& along,
                                                                                       TransformWith<TransformArgs...>&&
                                                                                           transArgs) const {
      const auto coeffsIndices = extractCoeffIndices(std::forward<Wrt<Args...>>(args));
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      maytransformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));

      return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(N, dNTransformed, std::get<0>(along.args),
                                                                           coeffsIndices);
    }

    /** \brief Function to forward the call of one spatial derivative in a single directions and two derivative wrt.
     * coefficients. You have to pass a along argument which specifies the direction where this derivative is applied */
    template <typename... Args, typename... TransformArgs, typename... AlongArgs,
              typename DomainTypeOrIntegrationPointIndex>
    requires(
        hasTwoCoeff<Wrt<Args...>>and DerivativeDirections::HasOneSpatialSingle<Wrt<Args...>, gridDim>and
            hasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex&
                                                                                                         localOrIpId,
                                                                                                     Wrt<Args...>&&
                                                                                                         args,
                                                                                                     Along<
                                                                                                         AlongArgs...>&&
                                                                                                         along,
                                                                                                     TransformWith<
                                                                                                         TransformArgs...>&&
                                                                                                         transArgs)
        const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      maytransformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));
      const auto coeffsIndices = extractCoeffIndices(std::forward<Wrt<Args...>>(args));

      const std::array<int, gridDim> counter
          = DerivativeDirections::countDynamicSpatialDerivativesInTuple<Wrt<Args...>, gridDim>(
              std::forward<Wrt<Args...>>(args));
      const int spatialIndex = DerivativeDirections::findSingleSpatial(counter);
      return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
          N, dNTransformed, std::get<0>(along.args), coeffsIndices, spatialIndex);
    }

    template <typename... TransformArgs>
    void maytransformDerivatives(const AnsatzFunctionJacobian& dNraw,
                                 TransformWith<TransformArgs...>&& transArgs) const {
      if constexpr (sizeof...(TransformArgs) > 0)
        dNTransformed = dNraw * std::get<0>(transArgs.args);
      else
        dNTransformed = dNraw;
    }

    template <typename DomainTypeOrIntegrationPointIndex>
    auto evaluateFunctionAndDerivativeWithIPorCoord(const DomainTypeOrIntegrationPointIndex& localOrIpId) const {
      if constexpr (std::is_same_v<DomainTypeOrIntegrationPointIndex, DomainType>) {
        AnsatzFunctionJacobian dN;
        impl().basis.evaluateJacobian(localOrIpId, dN);
        AnsatzFunctionType N;
        impl().basis.evaluateFunction(localOrIpId, N);
        return std::make_tuple(N, dN);
      } else if constexpr (std::numeric_limits<DomainTypeOrIntegrationPointIndex>::is_integer) {
        const AnsatzFunctionJacobian& dN = impl().basis.evaluateJacobian(localOrIpId);
        const AnsatzFunctionType& N      = impl().basis.evaluateFunction(localOrIpId);
        return std::make_tuple(std::ref(N), std::ref(dN));
      } else
        static_assert(std::is_same_v<DomainTypeOrIntegrationPointIndex,
                                     DomainType> or std::is_same_v<DomainTypeOrIntegrationPointIndex, int>,
                      "The argument you passed should be an id for the integration point or the point where the "
                      "derivative should be evaluated");
    }

    /** \brief Function to forward the call of one spatial derivative in all directions and no derivative wrt.
     * coefficients  */
    template <typename... Args, typename... TransformArgs, typename DomainTypeOrIntegrationPointIndex>
    requires(hasNoCoeff<Wrt<Args...>>and DerivativeDirections::HasOneSpatialAll<Wrt<Args...>, gridDim>and
                 hasevaluateDerivativeWRTSpaceAllImpl) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex&
                                                                                   localOrIpId,
                                                                               Wrt<Args...>&& args,
                                                                               TransformWith<TransformArgs...>&&
                                                                                   transArgs) const {
      const auto& [N, dN] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);

      maytransformDerivatives(dN, std::forward<TransformWith<TransformArgs...>>(transArgs));

      return impl().evaluateDerivativeWRTSpaceAllImpl(N, dNTransformed);
    }

    /** \brief Function to forward the call of one spatial derivative in a single directions and no derivative wrt.
     * coefficients  */
    template <typename... Args, typename... TransformArgs, typename DomainTypeOrIntegrationPointIndex>
    requires(
        hasNoCoeff<Wrt<Args...>>and DerivativeDirections::HasOneSpatialSingle<Wrt<Args...>, gridDim>and
            hasevaluateDerivativeWRTSpaceSingleImpl) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex&
                                                                                 localOrIpId,
                                                                             Wrt<Args...>&& args,
                                                                             TransformWith<TransformArgs...>&&
                                                                                 transArgs) const {
      const auto& [N, dN] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);

      maytransformDerivatives(dN, std::forward<TransformWith<TransformArgs...>>(transArgs));

      const std::array<int, gridDim> counter
          = DerivativeDirections::countDynamicSpatialDerivativesInTuple<Wrt<Args...>, gridDim>(
              std::forward<Wrt<Args...>>(args));
      auto singleSpatial = DerivativeDirections::findSingleSpatial<gridDim>(counter);
      return impl().evaluateDerivativeWRTSpaceSingleImpl(N, dNTransformed, singleSpatial);
    }

    /** \brief Function to forward the call of one spatial derivative and no derivative wrt. coefficients and No
     * transformation information  */
    template <typename... Args, typename... TransformArgs, typename DomainTypeOrIntegrationPointIndex>
    requires(hasNoCoeff<Wrt<Args...>>and DerivativeDirections::HasOneSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args) const {
      return evaluateDerivative(localOrIpId, std::forward<Wrt<Args...>>(args), transformWith());
    }

    /** \brief Function to forward the call of no spatial derivative and one derivative wrt. coefficients  */
    template <typename... Args, typename... Indices, typename DomainTypeOrIntegrationPointIndex>
    requires(hasSingleCoeff<Wrt<Args...>>and DerivativeDirections::HasNoSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args) const {
      const auto& [N, dN] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      const auto coeffsIndices = extractCoeffIndices(std::forward<Wrt<Args...>>(args));
      return impl().evaluateDerivativeWRTCoeffsImpl(N, dN, coeffsIndices);
    }

    /** \brief Function to forward the call of one spatial derivative in all directions and one derivative wrt.
     * coefficients  */
    template <typename... Args, typename... TransformArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasSingleCoeff<Wrt<Args...>>and DerivativeDirections::HasOneSpatialAll<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             TransformWith<TransformArgs...>&& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      maytransformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));
      const auto coeffsIndices = extractCoeffIndices(std::forward<Wrt<Args...>>(args));
      return impl().evaluateDerivativeWRTCoeffsANDSpatialImpl(N, dNTransformed, coeffsIndices);
    }

    /** \brief Function to forward the call of one spatial derivative in a single direction and one derivative wrt.
     * coefficients  */
    template <typename... Args, typename... TransformArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasSingleCoeff<Wrt<Args...>>and DerivativeDirections::HasOneSpatialSingle<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             TransformWith<TransformArgs...>&& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      maytransformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));

      const std::array<int, gridDim> counter
          = DerivativeDirections::countDynamicSpatialDerivativesInTuple<Wrt<Args...>, gridDim>(
              std::forward<Wrt<Args...>>(args));
      const int spatialIndex = DerivativeDirections::findSingleSpatial(counter);
      const auto coeffsIndices = extractCoeffIndices(std::forward<Wrt<Args...>>(args));
      return impl().evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(N, dNTransformed, coeffsIndices,
                                                                    spatialIndex);
    }

    /** \brief Function to forward the call of one spatial derivative and one derivative wrt. coefficients  and no
     * transformation arg*/
    template <typename... Args, typename... TransformArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasSingleCoeff<Wrt<Args...>>and DerivativeDirections::HasOneSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args) const {
      return evaluateDerivative(localOrIpId, std::forward<Wrt<Args...>>(args), transformWith());
    }

    auto viewOverIntegrationPoints() { return impl().basis.viewOverIntegrationPoints(); }

  private:
    auto tryCallSecondDerivativeWRTCoeffs(const auto& N, const auto& dN, const auto& along, const auto& coeffs) const {
      if constexpr (requires { impl().evaluateSecondDerivativeWRTCoeffs(N, dN, along, coeffs); })
        return impl().evaluateSecondDerivativeWRTCoeffs(N, dN, along, coeffs);
      else
        static_assert(
            requires { impl().evaluateSecondDerivativeWRTCoeffs(N, dN, along, coeffs); },
            "Your function does not have evaluateSecondDerivativeWRTCoeffs. Maybe this is on purpose and this would "
            "yield a zeroMatrix?");
    }

    mutable AnsatzFunctionJacobian dNTransformed;

    LocalFunctionImpl const& impl() const  // CRTP
    {
      return static_cast<LocalFunctionImpl const&>(*this);
    }
  };
}  // namespace Ikarus
