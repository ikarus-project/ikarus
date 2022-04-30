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

  template <typename LocalFunctionImpl>
  class LocalFunctionInterface;

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
                    "This type should not be instantiated. Check that your arguments satisfies the template below");
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
      const auto coeffIndicesOfArgs = Ikarus::DerivativeDirections::extractCoeffIndices(args);
      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<coeffIndicesOfArgs.size()>{}),
                            [&](auto&& i) { coeffsIndices[i]._data = coeffIndicesOfArgs[i]._data; });

      spatialPartialIndices = Ikarus::DerivativeDirections::extractSpatialPartialIndices(args);
    }

    // Constructor that does not calculate extractCoeffIndices and extractSpatialPartialIndices
    LocalFunctionEvaluationArgs(const DomainTypeOrIntegrationPointIndex& localOrIpId, const Wrt<WrtArgs...>& args,
                                const Along<AlongArgs...>& along, const TransformWith<TransformArgs...>& transArgs,
                                bool)
        : integrationPointOrIndex{localOrIpId}, wrtArgs{args}, alongArgs{along}, transformWithArgs{transArgs} {}

  public:
    auto extractSpatialOrFirstWrtArg() const {
      if constexpr (hasOneSpatial) {
        if constexpr (hasOneSpatialSingle)
          return extractWrtArgsWithGivenType<DerivativeDirections::SpatialPartial>();
        else if constexpr (hasOneSpatialAll)
          return extractWrtArgsWithGivenType<DerivativeDirections::SpatialAll>();
      } else
        return extractWrtArgs<0>();
    }

    auto extractSecondWrtArgOrFirstNonSpatial() const {
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

    template <typename DerivativeDirection>
    auto extractWrtArgsWithGivenType() const {
      if constexpr (std::is_same_v<DerivativeDirection, DerivativeDirections::SpatialPartial>) {
        auto wrtArg = wrt(DerivativeDirections::spatial(spatialPartialIndices));
        return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                           DomainTypeOrIntegrationPointIndex>(integrationPointOrIndex, wrtArg,
                                                                              alongArgs, transformWithArgs);
      } else if constexpr (std::is_same_v<DerivativeDirection, DerivativeDirections::SpatialAll>) {
        auto wrtArg = wrt(DerivativeDirections::spatialAll);
        return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                           DomainTypeOrIntegrationPointIndex>(integrationPointOrIndex, wrtArg,
                                                                              alongArgs, transformWithArgs);
      }
    }

    template <template <auto...> class DerivativeDirection>
    auto extractWrtArgsWithGivenType() const {
      using namespace Dune::Indices;
      if constexpr (Std::isTemplateSame_v<DerivativeDirections::TwoCoeff, DerivativeDirection>) {
        auto wrtArg = wrt(DerivativeDirections::coeff(coeffsIndices[_0][_0], coeffsIndices[_0][1],
                                                      coeffsIndices[_1][_0], coeffsIndices[_1][1]));
        return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                           DomainTypeOrIntegrationPointIndex>(integrationPointOrIndex, wrtArg,
                                                                              alongArgs, transformWithArgs);
      } else if constexpr (Std::isTemplateSame_v<DerivativeDirections::SingleCoeff, DerivativeDirection>) {
        auto wrtArg = wrt(DerivativeDirections::coeff(coeffsIndices[_0][_0], coeffsIndices[_0][1]));
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

  template <typename... WrtArgs, typename... OtherWrtArgs, typename... TransformArgs, typename... AlongArgs,
            typename DomainTypeOrIntegrationPointIndex>
  auto joinWRTArgs(
      const LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                        DomainTypeOrIntegrationPointIndex>& a,
      const LocalFunctionEvaluationArgs<Wrt<OtherWrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                        DomainTypeOrIntegrationPointIndex>& b) {
    auto wrt_lambda = [](auto... args) { return wrt(args...); };
    auto wrtArg     = std::apply(wrt_lambda, std::tuple_cat(a.wrtArgs.args, b.wrtArgs.args));
    return LocalFunctionEvaluationArgs<decltype(wrtArg), Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                       DomainTypeOrIntegrationPointIndex>(a.integrationPointOrIndex, wrtArg,
                                                                          a.alongArgs, a.transformWithArgs);
  }

  template <typename... WrtArgs, typename... OtherWrtArgs, typename... TransformArgs, typename... AlongArgs,
            typename DomainTypeOrIntegrationPointIndex>
  auto extractWrtArgsTwoCoeffsToSingleCoeff(
      const LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                        DomainTypeOrIntegrationPointIndex>& a) {
    using namespace Dune::Indices;
    const auto coeffArg = Std::getSpecialization<DerivativeDirections::TwoCoeff>(a.wrtArgs.args);
    auto wrtArg0        = wrt(DerivativeDirections::coeff(a.coeffsIndices[_0][_0], a.coeffsIndices[_0][1]));
    auto wrtArg1        = wrt(DerivativeDirections::coeff(a.coeffsIndices[_1][_0], a.coeffsIndices[_1][1]));

    return std::make_pair(
        LocalFunctionEvaluationArgs(a.integrationPointOrIndex, wrtArg0, a.alongArgs, a.transformWithArgs),
        LocalFunctionEvaluationArgs(a.integrationPointOrIndex, wrtArg1, a.alongArgs, a.transformWithArgs));
  }

  // This function returns the first two args and returns the spatial derivative argument always as first
  template <typename... WrtArgs, typename... OtherWrtArgs, typename... TransformArgs, typename... AlongArgs,
            typename DomainTypeOrIntegrationPointIndex>
  auto extractFirstTwoArgs(
      const LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                        DomainTypeOrIntegrationPointIndex>& a) {
    if constexpr (std::tuple_size_v<decltype(a.wrtArgs.args)> == 2) {
      if constexpr (DerivativeDirections::isSpatial<std::tuple_element_t<0, decltype(a.wrtArgs.args)>>)
        return std::make_pair(a.template extractWrtArgs<0>(), a.template extractWrtArgs<1>());
      else if constexpr (DerivativeDirections::isSpatial<std::tuple_element_t<1, decltype(a.wrtArgs.args)>>)
        return std::make_pair(a.template extractWrtArgs<1>(), a.template extractWrtArgs<0>());
    } else
      return extractWrtArgsTwoCoeffsToSingleCoeff(a);
  }

  template <typename... WrtArgs, typename... TransformArgs, typename... AlongArgs, typename... AlongArgsOther,
            typename DomainTypeOrIntegrationPointIndex>
  auto addAlong(const LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgs...>, TransformWith<TransformArgs...>,
                                                  DomainTypeOrIntegrationPointIndex>& args,
                const Along<AlongArgsOther...>& alongArgs) {
    auto newArgs = LocalFunctionEvaluationArgs<Wrt<WrtArgs...>, Along<AlongArgsOther...>,
                                               TransformWith<TransformArgs...>, DomainTypeOrIntegrationPointIndex>(
        args.integrationPointOrIndex, args.wrtArgs, alongArgs, args.transformWithArgs, false);

    using namespace Dune::Indices;
    std::get<1>(newArgs.coeffsIndices[_0]._data) = std::get<1>(args.coeffsIndices[_0]._data);
    std::get<1>(newArgs.coeffsIndices[_1]._data) = std::get<1>(args.coeffsIndices[_1]._data);
    newArgs.spatialPartialIndices                = args.spatialPartialIndices;

    return newArgs;
  }

  namespace Impl {

    template <typename LocalFunctionImpl_>
      requires(!std::is_arithmetic_v<LocalFunctionImpl_>)
    consteval int countUniqueNonArithmeticLeafNodesImpl() {


      if constexpr(Std::isSpecialization<std::tuple,typename LocalFunctionImpl_::Ids>::value) {
        constexpr auto predicate = []<typename Type>(Type ){return Type::value!=Ikarus::arithmetic;};
        return std::tuple_size_v<decltype(Std::unique(Std::filter(typename LocalFunctionImpl_::Ids(), predicate)))>;
      }
      else
        return 1;
    }


  template <typename LF> requires LocalFunction<LF>
   auto collectNonArithmeticLeafNodesImpl(const LF& a) {

//    static_assert(LocalFunction<LF>,"Only passing LocalFunctions allowed");
    if constexpr(IsBinaryExpr<LF>)
    return std::tuple_cat(collectNonArithmeticLeafNodesImpl(a.l()),collectNonArithmeticLeafNodesImpl(a.r()));
    else if constexpr(IsUnaryExpr<LF>)
      return std::make_tuple(collectNonArithmeticLeafNodesImpl(a.m()));
    else  if constexpr(IsArithmeticExpr<LF>)
    return std::make_tuple();
    else  if constexpr( IsNonArithmeticLeafNode<LF>)
      return std::make_tuple(std::cref(a));
    else
      static_assert("There are currently no other expressions. Thus you should not end up here.");


  }

  }  // namespace Impl

  template <typename LocalFunctionImpl_>
  consteval int countUniqueNonArithmeticLeafNodes(const LocalFunctionInterface<LocalFunctionImpl_>& a) {
    return Impl::countUniqueNonArithmeticLeafNodesImpl<LocalFunctionImpl_>();
  }


template <typename LF> requires LocalFunction<LF>
auto collectNonArithmeticLeafNodes(const LocalFunctionInterface<LF>& a) {

  return Std::makeNestedTupleFlatAndStoreReferences(Impl::collectNonArithmeticLeafNodesImpl(a.impl()));

}

template <typename LF> requires LocalFunction<LF>
 struct LocalFunctionLeafNodeCollection
 {
   LocalFunctionLeafNodeCollection(const LF& lf): leafNodes{collectNonArithmeticLeafNodes(lf)} {}

   template<std::size_t I>
   auto& coefficientsRef(Dune::index_constant<I> =Dune::index_constant<0UL>()) { return std::get<I>(leafNodes).coefficientsRef(); }
   template<std::size_t I>
   auto& basis(Dune::index_constant<I> =Dune::index_constant<0UL>()) { return std::get<I>(leafNodes).basis(); }

 private:
   decltype(collectNonArithmeticLeafNodes(std::declval<const LF&>())) leafNodes;
 };

template <typename LF> requires LocalFunction<LF>
    auto collectLeafNodeLocalFunctions(const LF& lf)
{
      return LocalFunctionLeafNodeCollection(lf);
}


  template <typename LocalFunctionImpl>
  class LocalFunctionInterface {
  public:
    using Traits     = LocalFunctionTraits<LocalFunctionImpl>;
    using DomainType = typename Traits::DomainType;
    static constexpr int gridDim =  Traits::gridDim;

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

  protected:
    /* Default implementation returns Zero expression if they are not overloaded */
    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateSecondDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const std::array<size_t, 2>& coeffsIndex,
                                               const Along<AlongArgs...>& alongArgs,
                                               const TransformWith<TransformArgs...>& transArgs) const {
      return Eigen::Matrix<typename LocalFunctionImpl::ctype, LocalFunctionImpl::correctionSize,
                           LocalFunctionImpl::correctionSize>::Zero();
    }

    /* Default implementation returns Zero expression if they are not overloaded */
    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const int spatialIndex, const Along<AlongArgs...>& alongArgs,
        const TransformWith<TransformArgs...>& transArgs) const {
      return Eigen::Matrix<typename LocalFunctionImpl::ctype, LocalFunctionImpl::correctionSize,
                           LocalFunctionImpl::correctionSize>::Zero();
    }

    /* Default implementation returns Zero expression if they are not overloaded */
    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    auto evaluateDerivativeWRTSpaceAllImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const TransformWith<TransformArgs...>& transArgs) const {
      return typename LocalFunctionImpl::Jacobian::Zero();
    }

    /* Default implementation returns Zero expression if they are not overloaded */
    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    auto evaluateDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           int coeffsIndex,
                                                           const TransformWith<TransformArgs...>& transArgs) const {
      return Eigen::Matrix<typename LocalFunctionImpl::ctype, LocalFunctionImpl::valueSize,
                           LocalFunctionImpl::correctionSize>::Zero();
    }

    /* Default implementation returns Zero expression if they are not overloaded  */
    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    auto evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const TransformWith<TransformArgs...>& transArgs) const {
      return std::array<Ikarus::DerivativeDirections::DerivativeNoOp, gridDim>();
    }

    /* Default implementation returns Zero expression if they are not overloaded  */
template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
auto evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
    const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
    const TransformWith<TransformArgs...>& transArgs) const {
  return Eigen::Matrix<typename LocalFunctionImpl::ctype, LocalFunctionImpl::valueSize,
                       LocalFunctionImpl::correctionSize>::Zero();
}
    /* Default implementation returns Zero expression if they are not overloaded  */
    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const Along<AlongArgs...>& alongArgs, const TransformWith<TransformArgs...>& transArgs) const {
      return std::array<Ikarus::DerivativeDirections::DerivativeNoOp, gridDim>();
    }

    /* Default implementation returns Zero expression if they are not overloaded  */
    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    auto evaluateDerivativeWRTSpaceSingleImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                         int spaceIndex,
                                                         const TransformWith<TransformArgs...>& transArgs) const {
      return typename Eigen::internal::plain_col_type<typename LocalFunctionImpl::Jacobian>::type::Zero();
    }



   private:
    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateDerivativeImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                       const LocalFunctionEvaluationArgs_& localFunctionArgs);

    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateFunctionImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                     const LocalFunctionEvaluationArgs_& localFunctionArgs);

    template <typename LF> requires LocalFunction<LF>
    friend auto collectNonArithmeticLeafNodes(const LocalFunctionInterface<LF>& a);

    constexpr LocalFunctionImpl const& impl() const  // CRTP
    {
      return static_cast<LocalFunctionImpl const&>(*this);
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
    using namespace Dune::Indices;
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
        if constexpr (decltype(localFunctionArgs.coeffsIndices[_0][_0])::value != LocalFunctionImpl::Ids::value)
          return DerivativeDirections::DerivativeNoOp();
        else if constexpr (localFunctionArgs.hasNoSpatial) {
          return f.impl().evaluateDerivativeWRTCoeffsImpl(localFunctionArgs.integrationPointOrIndex,
                                                          localFunctionArgs.coeffsIndices[_0][1],
                                                          localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialSingle) {
          return f.impl().evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
              localFunctionArgs.integrationPointOrIndex, localFunctionArgs.coeffsIndices[_0][1],
              localFunctionArgs.spatialPartialIndices, localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialAll) {
          return f.impl().evaluateDerivativeWRTCoeffsANDSpatialImpl(localFunctionArgs.integrationPointOrIndex,
                                                                    localFunctionArgs.coeffsIndices[_0][1],
                                                                    localFunctionArgs.transformWithArgs);
        }
      } else if constexpr (localFunctionArgs.hasTwoCoeff) {
        if constexpr (localFunctionArgs.hasNoSpatial) {
          return f.impl().evaluateSecondDerivativeWRTCoeffsImpl(
              localFunctionArgs.integrationPointOrIndex,
              {localFunctionArgs.coeffsIndices[_0][1], localFunctionArgs.coeffsIndices[_1][1]},
              localFunctionArgs.alongArgs, localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialSingle) {
          return f.impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
              localFunctionArgs.integrationPointOrIndex,
              {localFunctionArgs.coeffsIndices[_0][1], localFunctionArgs.coeffsIndices[_1][1]},
              localFunctionArgs.spatialPartialIndices, localFunctionArgs.alongArgs,
              localFunctionArgs.transformWithArgs);
        } else if constexpr (localFunctionArgs.hasOneSpatialAll) {
          return f.impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
              localFunctionArgs.integrationPointOrIndex,
              {localFunctionArgs.coeffsIndices[_0][1], localFunctionArgs.coeffsIndices[_1][1]},
              localFunctionArgs.alongArgs, localFunctionArgs.transformWithArgs);
        }
      }
    } else {
      return f.impl().template evaluateDerivativeOfExpression<localFunctionArgs.derivativeOrder>(localFunctionArgs);
    }
  }

  template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl>
  auto evaluateFirstOrderDerivativesImpl(const LocalFunctionInterface<LocalFunctionImpl>& f,
                                         const LocalFunctionEvaluationArgs_& localFunctionArgs) {
    if constexpr (localFunctionArgs.derivativeOrder == 3) {
      const auto argsForDx              = localFunctionArgs.extractSpatialOrFirstWrtArg();
      const auto [argsForDy, argsForDz] = extractWrtArgsTwoCoeffsToSingleCoeff(localFunctionArgs);
      auto dfdx                         = evaluateDerivativeImpl(f, argsForDx);
      auto dfdy                         = evaluateDerivativeImpl(f, argsForDy);
      auto dfdz                         = evaluateDerivativeImpl(f, argsForDz);
      return std::make_tuple(dfdx, dfdy, dfdz);
    } else if constexpr (localFunctionArgs.derivativeOrder == 2) {
      const auto [argsForDx, argsForDy] = extractFirstTwoArgs(localFunctionArgs);
      auto dfdx                         = evaluateDerivativeImpl(f, argsForDx);
      auto dfdy                         = evaluateDerivativeImpl(f, argsForDy);
      return std::make_tuple(dfdx, dfdy);
    }
  }

  template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl>
  auto evaluateSecondOrderDerivativesImpl(const LocalFunctionInterface<LocalFunctionImpl>& f,
                                          const LocalFunctionEvaluationArgs_& localFunctionArgs) {
    if constexpr (localFunctionArgs.derivativeOrder == 3) {
      const auto argsForDx              = localFunctionArgs.extractSpatialOrFirstWrtArg();
      const auto [argsForDy, argsForDz] = extractWrtArgsTwoCoeffsToSingleCoeff(localFunctionArgs);
      const auto argsForDxy             = joinWRTArgs(argsForDx, argsForDy);
      const auto argsForDxz             = joinWRTArgs(argsForDx, argsForDz);
      const auto df_dxy                 = evaluateDerivativeImpl(f, argsForDxy);
      const auto df_dxz                 = evaluateDerivativeImpl(f, argsForDxz);
      return std::make_tuple(df_dxy, df_dxz);
    } else
      static_assert(localFunctionArgs.derivativeOrder == 3);
  }

}  // namespace Ikarus
