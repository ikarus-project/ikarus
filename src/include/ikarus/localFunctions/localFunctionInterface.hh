//
// Created by alex on 3/17/22.
//

#pragma once
#include <concepts>

#include <ikarus/localBasis/localBasis.hh>
namespace Ikarus {
  namespace DerivativeDirections {
    constexpr std::integral_constant<int, -1> coeffs     = {};
    constexpr std::integral_constant<int, -2> spatialall = {};

    struct spatialPartial {
      using value = std::integral_constant<int, -3>;
      int index{};
    };

    spatialPartial spatial(int i);

    template <int Dim>
    struct Counter {
      int coeffDerivatives{};
      int spatialall{};
      int dynamicspatial{};
    };

    template <typename WrtType, int gridDim>
    consteval Counter<gridDim> countInTuple() {
      Counter<gridDim> counter{};
      Dune::Hybrid::forEach(
          Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(WrtType::args)>>()), [&](auto i) {
            using currentDerivType = std::tuple_element_t<i, typename WrtType::Args>;
            if constexpr (!std::is_same_v<currentDerivType, spatialPartial>) {
              constexpr int currentDerivValue = currentDerivType::value;
              constexpr bool isSpatialall     = spatialall.value == currentDerivValue;

              if constexpr (DerivativeDirections::coeffs.value == currentDerivValue)
                ++counter.coeffDerivatives;
              else if constexpr (isSpatialall)
                ++counter.spatialall;
            } else
              ++counter.dynamicspatial;
          });
      return counter;
    }

    template <typename WrtType, int gridDim>
    std::array<int, gridDim> countDynamicSpatialDerivativesInTuple(const WrtType& wrt) {
      std::array<int, gridDim> counter{};
      Dune::Hybrid::forEach(
          Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(WrtType::args)>>()), [&](auto i) {
            using currentDerivType = std::tuple_element_t<i, typename WrtType::Args>;
            if constexpr (std::is_same_v<currentDerivType, spatialPartial>) ++counter[std::get<i>(wrt.args).index];
          });
      return counter;
    }

    template <int gridDim>
    consteval bool hasNoSpatial(const std::array<int, gridDim>& spatialDerivs) {
      return std::ranges::count(spatialDerivs, 0) == gridDim;
    }

    template <int gridDim>
    consteval bool hasOneSpatial(const std::array<int, gridDim>& spatialDerivs) {
      return std::ranges::any_of(spatialDerivs, [](auto&& el) { return el != 0; }) > 0;
    }

    template <size_t gridDim>
    int findSingleSpatial(const std::array<int, gridDim>& spatialDerivs) {
      const auto element = std::ranges::find_if_not(spatialDerivs, [](auto&& el) { return el == 0; });
      return std::distance(spatialDerivs.begin(), element);
    }

    template <typename WrtType, int gridDim, int I>
    concept HasCoeff = (countInTuple<WrtType, gridDim>().coeffDerivatives == I);

    template <typename WrtType, int gridDim>
    concept HasNoSpatial
        = (countInTuple<WrtType, gridDim>().dynamicspatial == 0 and countInTuple<WrtType, gridDim>().spatialall == 0);

    template <typename WrtType, int gridDim>
    concept HasOneSpatialAll = countInTuple<WrtType, gridDim>().spatialall == 1;

    template <typename WrtType, int gridDim>
    concept HasOneSpatialSingle = (countInTuple<WrtType, gridDim>().dynamicspatial == 1);

    template <typename WrtType, int gridDim>
    concept HasOneSpatial = HasOneSpatialSingle<WrtType, gridDim> or HasOneSpatialAll<WrtType, gridDim>;

  }  // namespace DerivativeDirections

  template <typename... Args>
  struct CoeffIndices {
    std::array<std::common_type_t<std::remove_cvref_t<Args>...>, sizeof...(Args)> args{};
  };

  template <typename... Args>
  auto coeffIndices(Args&&... args) {
    return CoeffIndices<Args&&...>({std::forward<Args>(args)...});
  }

  template <typename... Args_>
  struct Wrt {
    using Args = std::tuple<std::remove_cvref_t<Args_>...>;
    std::tuple<std::remove_cvref_t<Args_>...> args;
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
    static constexpr bool hasevaluateDerivativeWRTSpaceAllImpl = HasevaluateDerivativeWRTSpaceAllImpl<LocalFunctionImpl>;
    static constexpr bool hasevaluateDerivativeWRTSpaceSingleImpl = HasevaluateDerivativeWRTSpaceSingleImpl<LocalFunctionImpl>;

    static constexpr int gridDim = Traits::gridDim;
    using TransformMatrix        = Eigen::Matrix<double, gridDim, gridDim>;

    template <typename WrtType, int I>
    static constexpr bool hasCoeff = DerivativeDirections::HasCoeff<WrtType, gridDim, I>;

    /** \brief Forward the binding to the local basis */
    template <typename IntegrationRule, typename... Ints>
    requires std::conjunction_v<std::is_convertible<int, Ints>...>
    void bind(IntegrationRule&& p_rule, Ikarus::Derivatives<Ints...>&& ints) {
      impl().basis.bind(std::forward<IntegrationRule>(p_rule), std::forward<Ikarus::Derivatives<Ints...>>(ints));
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
    template <typename... Args, typename... TransformArgs, typename... AlongArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(
        hasCoeff<Wrt<Args...>, 2>and DerivativeDirections::HasNoSpatial<Wrt<Args...>, gridDim>and
            hasSecondDerivativeWRTCoeffs) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                                  Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                                  TransformWith<TransformArgs...>&& transArgs,
                                                                  CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);

      // Check if a matrix is given to transform derivatives. Otherwise we do nothing
      if constexpr (sizeof...(TransformArgs) > 0) {
        AnsatzFunctionJacobian dN = (dNraw * std::get<0>(transArgs.args)).eval();
        return impl().evaluateSecondDerivativeWRTCoeffs(N, dN, std::get<0>(along.args), coeffsIndices.args);
      } else
        return impl().evaluateSecondDerivativeWRTCoeffs(N, dNraw, std::get<0>(along.args), coeffsIndices.args);
    }

    /** \brief Function to forward the call of no spatial derivative and two derivative wrt. coefficients.
     * You have to pass a along argument which specifies the direction where this derivative is applied
     * Specialization when no transformWith is passed
     */
    template <typename... Args, typename... AlongArgs, typename... Indices, typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 2>and DerivativeDirections::HasNoSpatial<Wrt<Args...>, gridDim>and
                 HasevaluateSecondDerivativeWRTCoeffs<
                     LocalFunctionImpl>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                                 Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                                 CoeffIndices<Indices...>&& coeffsIndices) const {
      return evaluateDerivative(localOrIpId, std::forward<Wrt<Args...>>(args), std::forward<Along<AlongArgs...>>(along),
                                transformWith(), std::forward<CoeffIndices<Indices...>>(coeffsIndices));
    }

    /** \brief Function to forward the call of one spatial derivative in all directions and two derivative wrt.
     * coefficients. You have to pass a along argument which specifies the direction where this derivative is applied */
    template <typename... Args, typename... TransformArgs, typename... AlongArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 2>and DerivativeDirections::HasOneSpatialAll<
             Wrt<Args...>, gridDim> and hasThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl)
        auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                             TransformWith<TransformArgs...>&& transArgs,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      maytransformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));

      return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(N, dNTransformed, std::get<0>(along.args),
                                                                           coeffsIndices.args);
    }

    /** \brief Function to forward the call of one spatial derivative in a single directions and two derivative wrt.
     * coefficients. You have to pass a along argument which specifies the direction where this derivative is applied */
    template <typename... Args, typename... TransformArgs, typename... AlongArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 2> and DerivativeDirections::HasOneSpatialSingle<
             Wrt<Args...>, gridDim> and hasevaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl)
        auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                             TransformWith<TransformArgs...>&& transArgs,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      maytransformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));

      const std::array<int, gridDim> counter
          = DerivativeDirections::countDynamicSpatialDerivativesInTuple<Wrt<Args...>, gridDim>(
              std::forward<Wrt<Args...>>(args));
      const int spatialIndex = DerivativeDirections::findSingleSpatial(counter);
      return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
          N, dNTransformed, std::get<0>(along.args), coeffsIndices.args, spatialIndex);
    }

    template <typename... TransformArgs>
    void maytransformDerivatives(const AnsatzFunctionJacobian& dNraw, TransformWith<TransformArgs...>&& transArgs) const {
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

    /** \brief Function to forward the call of one spatial derivative in all directions and no derivative wrt. coefficients  */
    template <typename... Args, typename... TransformArgs, typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 0>and DerivativeDirections::HasOneSpatialAll<
             Wrt<Args...>, gridDim> and hasevaluateDerivativeWRTSpaceAllImpl) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             TransformWith<TransformArgs...>&& transArgs) const {
      const auto& [N, dN] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);

      constexpr DerivativeDirections::Counter<gridDim> counterConstExpr
          = DerivativeDirections::countInTuple<Wrt<Args...>, gridDim>();

      maytransformDerivatives(dN, std::forward<TransformWith<TransformArgs...>>(transArgs));

      return impl().evaluateDerivativeWRTSpaceAllImpl(N, dNTransformed);
    }

    /** \brief Function to forward the call of one spatial derivative in a single directions and no derivative wrt. coefficients  */
    template <typename... Args, typename... TransformArgs, typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 0>and DerivativeDirections::HasOneSpatialSingle<
             Wrt<Args...>, gridDim> and hasevaluateDerivativeWRTSpaceSingleImpl) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             TransformWith<TransformArgs...>&& transArgs) const {
      const auto& [N, dN] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);

      constexpr DerivativeDirections::Counter<gridDim> counterConstExpr
          = DerivativeDirections::countInTuple<Wrt<Args...>, gridDim>();

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
    requires(hasCoeff<Wrt<Args...>, 0>and DerivativeDirections::HasOneSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args) const {
      return evaluateDerivative(localOrIpId, std::forward<Wrt<Args...>>(args), transformWith());
    }

    /** \brief Function to forward the call of no spatial derivative and one derivative wrt. coefficients  */
    template <typename... Args, typename... Indices, typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 1>and DerivativeDirections::HasNoSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dN] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      return impl().evaluateDerivativeWRTCoeffsImpl(N, dN, coeffsIndices.args[0]);
    }

    /** \brief Function to forward the call of one spatial derivative in all directions and one derivative wrt. coefficients  */
    template <typename... Args, typename... TransformArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 1>and DerivativeDirections::HasOneSpatialAll<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             TransformWith<TransformArgs...>&& transArgs,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      maytransformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));

      return impl().evaluateDerivativeWRTCoeffsANDSpatialImpl(N, dNTransformed, coeffsIndices.args[0]);
    }

    /** \brief Function to forward the call of one spatial derivative in a single direction and one derivative wrt. coefficients  */
    template <typename... Args, typename... TransformArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 1>and DerivativeDirections::HasOneSpatialSingle<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             TransformWith<TransformArgs...>&& transArgs,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      maytransformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));

      const std::array<int, gridDim> counter
            = DerivativeDirections::countDynamicSpatialDerivativesInTuple<Wrt<Args...>, gridDim>(
                std::forward<Wrt<Args...>>(args));
        const int spatialIndex = DerivativeDirections::findSingleSpatial(counter);
        return impl().evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(N, dNTransformed, coeffsIndices.args[0],
                                                                      spatialIndex);
    }

    /** \brief Function to forward the call of one spatial derivative and one derivative wrt. coefficients  and no
     * transformation arg*/
    template <typename... Args, typename... TransformArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 1>and DerivativeDirections::HasOneSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      return evaluateDerivative(localOrIpId, std::forward<Wrt<Args...>>(args), transformWith(),
                                std::forward<CoeffIndices<Indices...>>(coeffsIndices));
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
