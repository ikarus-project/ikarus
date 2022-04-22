//
// Created by alex on 3/17/22.
//

#pragma once
#include <concepts>

#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

  template <typename... Args_>
  struct Wrt {
    using Args = std::tuple<std::remove_cvref_t<Args_>...>;
    Args args;
  };

  template <typename... Args>
  auto wrt(Args&&... args) {
    return Wrt<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
  }

  template <typename... Args_>
  struct Along {
    using Args = std::tuple<Args_...>;
    Args args;
  };

  template <typename... Args>
  auto along(Args&&... args) {
    return Along<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
  }

  template <typename... Args_>
  struct TransformWith {
    using Args = std::tuple<Args_...>;
    Args args;
  };

  template <typename... Args>
  auto transformWith(Args&&... args) {
    return TransformWith<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
  }

  namespace DerivativeDirections {
    static struct Spatialall { using value = std::integral_constant<int, -2>; } spatialall;

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
      std::array<size_t, 2> index{};
    };

    SpatialPartial spatial(size_t i);
    SingleCoeff coeff(size_t i);
    TwoCoeff coeff(size_t i, size_t j);

    template <int Dim>
    struct Counter {
      int coeffDerivatives{};
      int spatialall{};
      int dynamicspatial{};
    };

    template <int gridDim, int coeffDerivs = 2>
    struct NewCounter {
      std::array<int, gridDim> spatialall{};
      int dynamicspatial{};
    };

    template <typename WrtType>
    auto extractSpatialPartialIndices(WrtType&& wrt) {
      if constexpr (Std::hasType<SpatialPartial, typename std::remove_reference_t<WrtType>::Args>::value)
        return std::get<SpatialPartial>(wrt.args).index;  // returns single int
      else
        return std::array<int, 0>();  // signals no SpatialPartial derivative found
    }

    struct ConstExprCounter {
      int singleCoeffDerivs{};
      int twoCoeffDerivs{};
      int spatialDerivs{};
      int spatialAll{};

      consteval int orderOfDerivative() const {
        return singleCoeffDerivs + 2 * twoCoeffDerivs + spatialDerivs + spatialAll;
      }
    };

    template <typename WrtType>
    consteval ConstExprCounter countDerivativesType() {
      ConstExprCounter counter{};
      using Tuple               = typename WrtType::Args;
      counter.singleCoeffDerivs = Ikarus::Std::countType<Tuple, SingleCoeff>();
      counter.twoCoeffDerivs    = Ikarus::Std::countType<Tuple, TwoCoeff>();
      counter.spatialDerivs     = Ikarus::Std::countType<Tuple, SpatialPartial>();
      counter.spatialAll        = Ikarus::Std::countType<Tuple, Spatialall>();
      return counter;
    }

    template <typename WrtType>
    auto extractCoeffIndices(WrtType&& wrt) {
      if constexpr (Std::hasType<SingleCoeff, typename std::remove_reference_t<WrtType>::Args>::value)
        return std::get<SingleCoeff>(wrt.args).index;  // returns single int
      else if constexpr (Std::hasType<TwoCoeff, typename std::remove_reference_t<WrtType>::Args>::value)
        return std::get<TwoCoeff>(wrt.args).index;  // return std::array<size_t,2>
      else
        return std::array<int, 0>();  // signals no coeff derivative found
    }

    template <size_t gridDim>
    int findSingleSpatial(const std::array<int, gridDim>& spatialDerivs) {
      const auto element = std::ranges::find_if_not(spatialDerivs, [](auto&& el) { return el == 0; });
      return std::distance(spatialDerivs.begin(), element);
    }

    template <typename WrtType>
    concept HasTwoCoeff = (countDerivativesType<WrtType>().twoCoeffDerivs == 1);

    template <typename WrtType>
    concept HasSingleCoeff = (countDerivativesType<WrtType>().singleCoeffDerivs == 1);

    template <typename WrtType>
    concept HasNoCoeff = (countDerivativesType<WrtType>().singleCoeffDerivs == 0
                          and countDerivativesType<WrtType>().twoCoeffDerivs == 0);

    template <typename WrtType>
    concept HasNoSpatial
        = (countDerivativesType<WrtType>().spatialDerivs == 0 and countDerivativesType<WrtType>().spatialAll == 0);

    template <typename WrtType>
    concept HasOneSpatialAll = countDerivativesType<WrtType>()
    .spatialAll == 1;

    template <typename WrtType>
    concept HasOneSpatialSingle = countDerivativesType<WrtType>()
    .spatialDerivs == 1;

    template <typename WrtType>
    concept HasOneSpatial = HasOneSpatialSingle<WrtType> or HasOneSpatialAll<WrtType>;

  }  // namespace DerivativeDirections

  template <typename Derived>
  struct LocalFunctionTraits;

  template <typename LocalFunctionImpl>
  concept HasevaluateSecondDerivativeWRTCoeffs
      = requires(LocalFunctionImpl func) {
          func.evaluateSecondDerivativeWRTCoeffs(
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>(),
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AlongType>(),
              std::declval<std::array<size_t, 2>>());
        };

  template <typename LocalFunctionImpl>
  concept HasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl
      = requires(LocalFunctionImpl func) {
          func.evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>(),
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AlongType>(),
              std::declval<std::array<size_t, 2>>());
        };

  template <typename LocalFunctionImpl>
  concept HasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl
      = requires(LocalFunctionImpl func) {
          func.evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>(),
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AlongType>(),
              std::declval<std::array<size_t, 2>>(), std::declval<int>());
        };

  template <typename LocalFunctionImpl>
  concept HasevaluateDerivativeWRTSpaceAllImpl
      = requires(LocalFunctionImpl func) {
          func.evaluateDerivativeWRTSpaceAllImpl(
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionType>(),
              std::declval<typename LocalFunctionTraits<LocalFunctionImpl>::AnsatzFunctionJacobian>());
        };

  template <typename DomainTypeOrIntegrationPointIndex, typename DomainType>
  concept IsIntegrationPointIndexOrIntegrationPointPosition
      = std::is_same_v<DomainTypeOrIntegrationPointIndex,
                       DomainType> or std::numeric_limits<DomainTypeOrIntegrationPointIndex>::is_integer;

  template <typename TypeListOne, typename TypeListTwo, typename TypeListThree,
            typename DomainTypeOrIntegrationPointIndex>
  struct LocalFunctionEvaluationArgs {
  public:
    LocalFunctionEvaluationArgs(const DomainTypeOrIntegrationPointIndex&, [[maybe_unused]] const TypeListOne& l1,
                                [[maybe_unused]] const TypeListTwo& l2, [[maybe_unused]] const TypeListThree& l3) {
      static_assert(!sizeof(TypeListOne),
                    "This type should not be instantiated. check that your arguments satisfies the template below");
    }
  };

  template <typename... WrtArgs, typename... TransformArgs, typename... AlongArgs,
            typename DomainTypeOrIntegrationPointIndex>
  struct LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                     DomainTypeOrIntegrationPointIndex> {
    LocalFunctionEvaluationArgs(const DomainTypeOrIntegrationPointIndex& localOrIpId, const Wrt<WrtArgs...>& args,
                                const Along<AlongArgs...>& along, const TransformWith<TransformArgs...>& transArgs)
        : integrationPointOrIndex{localOrIpId}, wrtArgs{args}, alongArgs{along}, transformWithArgs{transArgs} {
      coeffsIndices         = Ikarus::DerivativeDirections::extractCoeffIndices(args);
      spatialPartialIndices = Ikarus::DerivativeDirections::extractSpatialPartialIndices(args);
    }

    static constexpr DerivativeDirections::ConstExprCounter derivativeCounter
        = DerivativeDirections::countDerivativesType<Wrt<WrtArgs...>>();
    static constexpr int derivativeOrder
        = DerivativeDirections::countDerivativesType<Wrt<WrtArgs...>>().orderOfDerivative();

    static constexpr bool hasTwoCoeff         = DerivativeDirections::HasTwoCoeff<Wrt<WrtArgs...>>;
    static constexpr bool hasSingleCoeff      = DerivativeDirections::HasSingleCoeff<Wrt<WrtArgs...>>;
    static constexpr bool hasNoCoeff          = DerivativeDirections::HasNoCoeff<Wrt<WrtArgs...>>;
    static constexpr bool hasNoSpatial        = DerivativeDirections::HasNoSpatial<Wrt<WrtArgs...>>;
    static constexpr bool hasOneSpatialAll    = DerivativeDirections::HasOneSpatialAll<Wrt<WrtArgs...>>;
    static constexpr bool hasOneSpatialSingle = DerivativeDirections::HasOneSpatialSingle<Wrt<WrtArgs...>>;
    static constexpr bool hasOneSpatial       = hasOneSpatialAll or hasOneSpatialSingle;

    DomainTypeOrIntegrationPointIndex integrationPointOrIndex{};

    const Wrt<WrtArgs&&...>& wrtArgs;
    const Along<AlongArgs&&...>& alongArgs;
    const TransformWith<TransformArgs&&...>& transformWithArgs;
    decltype(DerivativeDirections::extractCoeffIndices<Wrt<WrtArgs...>>(std::declval<Wrt<WrtArgs...>>())) coeffsIndices;
    decltype(DerivativeDirections::extractSpatialPartialIndices<Wrt<WrtArgs...>>(
        std::declval<Wrt<WrtArgs...>>())) spatialPartialIndices;
  };

  template <typename LocalFunctionImpl>
  class LocalFunctionInterface {
  public:
    using Traits                 = LocalFunctionTraits<LocalFunctionImpl>;
    using DomainType             = typename Traits::DomainType;

    // Check the capabilities of derived implementation
    static constexpr bool hasSecondDerivativeWRTCoeffs = HasevaluateSecondDerivativeWRTCoeffs<LocalFunctionImpl>;
    static constexpr bool hasThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl
        = HasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl<LocalFunctionImpl>;
    static constexpr bool hasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl
        = HasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl<LocalFunctionImpl>;
    static constexpr bool hasevaluateDerivativeWRTSpaceAllImpl
        = HasevaluateDerivativeWRTSpaceAllImpl<LocalFunctionImpl>;

    static constexpr int gridDim = Traits::gridDim;
    using TransformMatrix        = Eigen::Matrix<double, gridDim, gridDim>;

    template <typename WrtType>
    static constexpr bool hasTwoCoeff = DerivativeDirections::HasTwoCoeff<WrtType>;
    template <typename WrtType>
    static constexpr bool hasSingleCoeff = DerivativeDirections::HasSingleCoeff<WrtType>;
    template <typename WrtType>
    static constexpr bool hasNoCoeff = DerivativeDirections::HasNoCoeff<WrtType>;

    /** \brief Forward the binding to the local basis */
    template <typename IntegrationRule, typename... Ints>
      requires std::conjunction_v<std::is_convertible<int, Ints>
                                  ...> void
      bind(IntegrationRule&& p_rule, Impl::Derivatives<Ints...>&& ints) {
      impl().basis.bind(std::forward<IntegrationRule>(p_rule), std::forward<Impl::Derivatives<Ints...>>(ints));
    }

    /** \brief Return the function value*/
    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
      requires IsIntegrationPointIndexOrIntegrationPointPosition<DomainTypeOrIntegrationPointIndex, DomainType>
    auto evaluateFunction(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                          const TransformWith<TransformArgs...>& transArgs = transformWith()) const {
      const LocalFunctionEvaluationArgs evalArgs(ipIndexOrPosition, wrt(), along(), transArgs);
      return impl().evaluateFunctionImpl(evalArgs.integrationPointOrIndex, evalArgs.transformWithArgs);
    }

    template <typename... WrtArgs, typename... TransformArgs, typename... AlongArgs,
              typename DomainTypeOrIntegrationPointIndex>
    auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId, Wrt<WrtArgs...>&& args,
                            Along<AlongArgs...>&& along,
                            TransformWith<TransformArgs...>&& transArgs = transformWith()) const {
      const LocalFunctionEvaluationArgs evalArgs(localOrIpId, std::forward<Wrt<WrtArgs...>>(args),
                                                 std::forward<Along<AlongArgs...>>(along), transArgs);
      return evaluateDerivativeImpl(*this, evalArgs);
    }

    template <typename... WrtArgs, typename... TransformArgs, typename DomainTypeOrIntegrationPointIndex>
      requires(
          (hasNoCoeff<Wrt<
               WrtArgs...>> and (DerivativeDirections::HasOneSpatialSingle<Wrt<WrtArgs...>> or DerivativeDirections::HasOneSpatialAll<Wrt<WrtArgs...>>))
          or (hasSingleCoeff<Wrt<WrtArgs...>> and DerivativeDirections::HasNoSpatial<Wrt<WrtArgs...>>)
          or (hasSingleCoeff<Wrt<WrtArgs...>> and DerivativeDirections::HasOneSpatialSingle<Wrt<WrtArgs...>>)
          or (hasSingleCoeff<Wrt<WrtArgs...>> and DerivativeDirections::HasOneSpatialAll<Wrt<WrtArgs...>>)
          or (hasTwoCoeff<Wrt<
                  WrtArgs...>> and DerivativeDirections::HasOneSpatialAll<Wrt<WrtArgs...>> and hasThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl)
          or (hasTwoCoeff<Wrt<
                  WrtArgs...>> and DerivativeDirections::HasOneSpatialSingle<Wrt<WrtArgs...>> and hasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl))
    auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId, Wrt<WrtArgs...>&& args,
                            TransformWith<TransformArgs...>&& transArgs = transformWith()) const {
      return evaluateDerivative(localOrIpId, std::forward<Wrt<WrtArgs...>>(args), along(),
                                std::forward<TransformWith<TransformArgs...>>(transArgs));
    }

    auto viewOverIntegrationPoints() { return impl().basis().viewOverIntegrationPoints(); }

  private:
    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateDerivativeImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                       const LocalFunctionEvaluationArgs_& localFunctionArgs);

    LocalFunctionImpl const& impl() const  // CRTP
    {
      return static_cast<LocalFunctionImpl const&>(*this);
    }
  };

  template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl>
  auto evaluateDerivativeImpl(const LocalFunctionInterface<LocalFunctionImpl>& f,
                              const LocalFunctionEvaluationArgs_& localFunctionArgs) {
    if constexpr (LocalFunctionImpl::isLeaf) {
      if constexpr (localFunctionArgs.hasNoCoeff) {
        if constexpr (localFunctionArgs.hasOneSpatialSingle) {
          return f.impl().evaluateDerivativeWRTSpaceSingleImpl(localFunctionArgs.integrationPointOrIndex,
                                                               localFunctionArgs.spatialPartialIndices,
                                                               localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialAll) {
          return f.impl().evaluateDerivativeWRTSpaceAllImpl(localFunctionArgs.integrationPointOrIndex,
                                                            localFunctionArgs.transformWithArgs);
        }
      } else if constexpr (localFunctionArgs.hasSingleCoeff) {
        if constexpr (localFunctionArgs.hasNoSpatial) {
          return f.impl().evaluateDerivativeWRTCoeffsImpl(localFunctionArgs.integrationPointOrIndex,
                                                          localFunctionArgs.coeffsIndices,
                                                          localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialSingle) {
          return f.impl().evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
              localFunctionArgs.integrationPointOrIndex, localFunctionArgs.coeffsIndices,
              localFunctionArgs.spatialPartialIndices, localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialAll) {
          return f.impl().evaluateDerivativeWRTCoeffsANDSpatialImpl(localFunctionArgs.integrationPointOrIndex,
                                                                    localFunctionArgs.coeffsIndices,
                                                                    localFunctionArgs.transformWithArgs);
        }
      } else if constexpr (localFunctionArgs.hasTwoCoeff) {
        if constexpr (localFunctionArgs.hasNoSpatial) {
          return f.impl().evaluateSecondDerivativeWRTCoeffsImpl(
              localFunctionArgs.integrationPointOrIndex, localFunctionArgs.coeffsIndices, localFunctionArgs.alongArgs,
              localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialSingle) {
          return f.impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
              localFunctionArgs.integrationPointOrIndex, localFunctionArgs.coeffsIndices,
              localFunctionArgs.spatialPartialIndices, localFunctionArgs.alongArgs,
              localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialAll) {
          return f.impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
              localFunctionArgs.integrationPointOrIndex, localFunctionArgs.coeffsIndices, localFunctionArgs.alongArgs,
              localFunctionArgs.transformWithArgs);
        }
      }
    } else {
      return f.impl().template evaluateDerivativeOfExpression<localFunctionArgs.derivativeOrder>(localFunctionArgs);
    }
  }

}  // namespace Ikarus
