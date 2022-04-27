//
// Created by alex on 3/17/22.
//

#pragma once
#include "meta.hh"

#include <concepts>

#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

  template <typename Derived>
  struct LocalFunctionTraits;

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
    template <typename, typename, typename, typename>
    friend class LocalFunctionEvaluationArgs;

    LocalFunctionEvaluationArgs(const DomainTypeOrIntegrationPointIndex& localOrIpId, const Wrt<WrtArgs...>& args,
                                const Along<AlongArgs...>& along, const TransformWith<TransformArgs...>& transArgs)
        : integrationPointOrIndex{localOrIpId}, wrtArgs{args}, alongArgs{along}, transformWithArgs{transArgs} {
      coeffsIndices         = Ikarus::DerivativeDirections::extractCoeffIndices(args);
      spatialPartialIndices = Ikarus::DerivativeDirections::extractSpatialPartialIndices(args);
    }


    // Constructor that does not calculate extractCoeffIndices and extractSpatialPartialIndices
    LocalFunctionEvaluationArgs(const DomainTypeOrIntegrationPointIndex& localOrIpId, const Wrt<WrtArgs...>& args,
                                const Along<AlongArgs...>& along, const TransformWith<TransformArgs...>& transArgs,
                                bool)
        : integrationPointOrIndex{localOrIpId}, wrtArgs{args}, alongArgs{along}, transformWithArgs{transArgs} {}

  public:
    auto extractSpatialOrFirstWrtArg() const{
      if constexpr (hasOneSpatial) {
        if constexpr (hasOneSpatialSingle)
          return extractWrtArgsWithGivenType<DerivativeDirections::SpatialPartial>();
         else if constexpr (hasOneSpatialAll)
          return extractWrtArgsWithGivenType<DerivativeDirections::Spatialall>();
      } else
        return extractWrtArgs<0>();
    }

    auto extractSecondWrtArgOrFirstNonSpatial()const {
      if constexpr (!hasOneSpatial)
        return extractWrtArgs<0>();
      else {
        static_assert(!hasNoCoeff, " There needs to be a coeff wrt argument!");

        if constexpr (hasSingleCoeff) {
          return extractWrtArgsWithGivenType<DerivativeDirections::SingleCoeff>();
        } else if constexpr (hasTwoCoeff) {
          return extractWrtArgsWithGivenType<DerivativeDirections::TwoCoeff>();
        }
      }
    }

    template <std::size_t... I>
    auto extractWrtArgs() const {
      auto wrt_lambda = [](auto... args) { return wrt(args...); };
      auto wrtArg
          = std::apply(wrt_lambda, Std::makeTupleFromTupleIndices<const decltype(wrtArgs.args)&, I...>(wrtArgs.args));
      return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                         DomainTypeOrIntegrationPointIndex>(integrationPointOrIndex, wrtArg, alongArgs,
                                                                            transformWithArgs);
    }

    auto extractWrtArgsTwoCoeffsToSingleCoeff() const
    {

      const auto coeffArg = std::get<DerivativeDirections::TwoCoeff>(wrtArgs.args);
      auto wrtArg0            = wrt(DerivativeDirections::coeff(coeffsIndices[0]));
      auto wrtArg1            = wrt(DerivativeDirections::coeff(coeffsIndices[1]));
      using ReturnedArgs = LocalFunctionEvaluationArgs<decltype(wrtArg0), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                                       DomainTypeOrIntegrationPointIndex>;

      return std::make_pair(ReturnedArgs(integrationPointOrIndex, wrtArg0,   alongArgs, transformWithArgs),ReturnedArgs(integrationPointOrIndex, wrtArg1,    alongArgs, transformWithArgs));
    }

    template<typename DerivativeDirection>
    auto extractWrtArgsWithGivenType()const
    {
      if constexpr(std::is_same_v<DerivativeDirection,DerivativeDirections::TwoCoeff>) {
        auto wrtArg = wrt(DerivativeDirections::coeff(coeffsIndices[0], coeffsIndices[1]));
        return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                           DomainTypeOrIntegrationPointIndex>(integrationPointOrIndex, wrtArg,
                                                                              alongArgs, transformWithArgs);
      }else if constexpr(std::is_same_v<DerivativeDirection,DerivativeDirections::SingleCoeff>) {
        auto wrtArg = wrt(DerivativeDirections::coeff(coeffsIndices));
        return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                           DomainTypeOrIntegrationPointIndex>(integrationPointOrIndex, wrtArg,
                                                                              alongArgs, transformWithArgs);
      }else if constexpr(std::is_same_v<DerivativeDirection,DerivativeDirections::SpatialPartial>) {
        auto wrtArg = wrt(DerivativeDirections::spatial(spatialPartialIndices));
        return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                           DomainTypeOrIntegrationPointIndex>(integrationPointOrIndex, wrtArg,
                                                                              alongArgs, transformWithArgs);
      }else if constexpr(std::is_same_v<DerivativeDirection,DerivativeDirections::Spatialall>) {
        auto wrtArg = wrt(DerivativeDirections::spatialall);
        return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                           DomainTypeOrIntegrationPointIndex>(integrationPointOrIndex, wrtArg,
                                                                              alongArgs, transformWithArgs);
      }
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

    const Wrt<WrtArgs&&...> wrtArgs;
    const Along<AlongArgs&&...> alongArgs;
    const TransformWith<TransformArgs&&...> transformWithArgs;
    decltype(DerivativeDirections::extractCoeffIndices<Wrt<WrtArgs...>>(std::declval<Wrt<WrtArgs...>>())) coeffsIndices;
    decltype(DerivativeDirections::extractSpatialPartialIndices<Wrt<WrtArgs...>>(
        std::declval<Wrt<WrtArgs...>>())) spatialPartialIndices;
  };

template <typename... WrtArgs,typename... OtherWrtArgs, typename... TransformArgs, typename... AlongArgs,
    typename DomainTypeOrIntegrationPointIndex>
 auto joinWRTArgs(const LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                                          DomainTypeOrIntegrationPointIndex>& a,
                  const LocalFunctionEvaluationArgs<Wrt<OtherWrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                                          DomainTypeOrIntegrationPointIndex>& b)
{
  auto wrt_lambda = [](auto... args) { return wrt(args...); };
  auto wrtArg
      = std::apply(wrt_lambda, std::tuple_cat(a.wrtArgs.args,b.wrtArgs.args));
  return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                     DomainTypeOrIntegrationPointIndex>(a.integrationPointOrIndex, wrtArg, a.alongArgs,
                                                                        a.transformWithArgs);
}


template <typename... WrtArgs, typename... TransformArgs, typename... AlongArgs,typename... AlongArgsOther,
    typename DomainTypeOrIntegrationPointIndex>
 auto createWithAlong(const LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                                              DomainTypeOrIntegrationPointIndex>& args, const Along<AlongArgsOther...>& alongArgs) {
  auto newArgs = LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgsOther...>,
                                             TransformWith<TransformArgs...>, DomainTypeOrIntegrationPointIndex>(
      args.integrationPointOrIndex, args.wrtArgs, alongArgs, args.transformWithArgs, false);

  newArgs.coeffsIndices         = args.coeffsIndices;
  newArgs.spatialPartialIndices = args.spatialPartialIndices;

  return newArgs;
}

  template <typename LocalFunctionImpl>
  class LocalFunctionInterface {
  public:
    using Traits     = LocalFunctionTraits<LocalFunctionImpl>;
    using DomainType = typename Traits::DomainType;

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
      return evaluateFunctionImpl(*this, evalArgs);
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

    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateFunctionImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                     const LocalFunctionEvaluationArgs_& localFunctionArgs);

    LocalFunctionImpl const& impl() const  // CRTP
    {
      return static_cast<LocalFunctionImpl const&>(*this);
    }

   protected:
   /*
    * Default implementation returns Zero expression if they are not overloaded
    */
    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateSecondDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const std::array<size_t, 2>& coeffsIndex,
                                               const Along<AlongArgs...>& alongArgs,
                                               const TransformWith<TransformArgs...>& transArgs) const {
      return Eigen::Matrix<typename LocalFunctionImpl::ctype, LocalFunctionImpl::correctionSize,LocalFunctionImpl::correctionSize>::Zero();
    }

       /*
    * Default implementation returns Zero expression if they are not overloaded
    */
    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const int spatialIndex, const Along<AlongArgs...>& alongArgs,
        const TransformWith<TransformArgs...>& transArgs) const {
      return Eigen::Matrix<typename LocalFunctionImpl::ctype, LocalFunctionImpl::correctionSize,LocalFunctionImpl::correctionSize>::Zero();
    }

  };

  template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl>
  auto evaluateFunctionImpl(const LocalFunctionInterface<LocalFunctionImpl>& f,
                            const LocalFunctionEvaluationArgs_& localFunctionArgs) {
    if constexpr (LocalFunctionImpl::isLeaf)
      return f.impl().evaluateFunctionImpl(localFunctionArgs.integrationPointOrIndex,
                                           localFunctionArgs.transformWithArgs);
    else {
      return f.impl().evaluateValueOfExpression(localFunctionArgs);
    }
  }

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

template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl>
auto evaluateFirstOrderDerivativesImpl(const LocalFunctionInterface<LocalFunctionImpl>& f,
                            const LocalFunctionEvaluationArgs_& localFunctionArgs) {

  return std::make_tuple()
}

}  // namespace Ikarus
