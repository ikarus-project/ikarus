//
// Created by alex on 3/17/22.
//

#pragma once
#include <concepts>
namespace Ikarus {
  namespace DerivativeDirections {
    constexpr std::integral_constant<int, -1> coeffs = {};

    template <int i>
    requires(i >= 0 and i < 3) using spatial             = std::integral_constant<int, i>;
    constexpr spatial<0> spatial0                        = {};
    constexpr spatial<1> spatial1                        = {};
    constexpr spatial<2> spatial2                        = {};
    constexpr std::integral_constant<int, -2> spatialall = {};

    template <int Dim>
    struct Counter {
      int coeffDerivatives{};
      std::array<int, Dim> spatialDerivatives{};
      int spatialall{};
    };

    template <typename WrtType, int gridDim>
    consteval Counter<gridDim> countInTuple() {
      Counter<gridDim> counter{};
      Dune::Hybrid::forEach(
          Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(WrtType::args)>>()), [&](auto i) {
            constexpr int currentDerivType = std::tuple_element_t<i, typename WrtType::Args>::value;
            constexpr bool isSpatial       = spatial0.value == currentDerivType or spatial1.value == currentDerivType
                                       or spatial2.value == currentDerivType or spatialall.value == currentDerivType;

            if constexpr (DerivativeDirections::coeffs.value == currentDerivType)
              ++counter.coeffDerivatives;
            else if constexpr (isSpatial) {
              if constexpr (DerivativeDirections::spatialall.value != currentDerivType)
                ++counter.spatialDerivatives[currentDerivType];
              else
                ++counter.spatialall;
            } else
              static_assert((DerivativeDirections::coeffs.value == currentDerivType) or isSpatial,
                            "You are only allowed to pass coeffs and spatial<n> to wrt");
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
    consteval int findSingleSpatial(const std::array<int, gridDim>& spatialDerivs) {
      const auto element = std::ranges::find_if_not(spatialDerivs, [](auto&& el) { return el == 0; });
      return std::distance(spatialDerivs.begin(), element);
    }

    template <typename WrtType, int gridDim, int I>
    concept HasCoeff = (countInTuple<WrtType, gridDim>().coeffDerivatives == I);

    template <typename WrtType, int gridDim>
    concept HasNoSpatial = (hasNoSpatial<gridDim>(countInTuple<WrtType, gridDim>().spatialDerivatives)
                            and countInTuple<WrtType, gridDim>().spatialall == 0);

    template <typename WrtType, int gridDim>
    concept HasOneSpatial = (hasOneSpatial<gridDim>(countInTuple<WrtType, gridDim>().spatialDerivatives)
                             or countInTuple<WrtType, gridDim>().spatialall > 0);

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
  class LocalFunctionInterface {
  public:
    using Traits                 = LocalFunctionTraits<LocalFunctionImpl>;
    using DomainType             = typename Traits::DomainType;
    using FunctionReturnType     = typename Traits::FunctionReturnType;
    using AnsatzFunctionType     = typename Traits::AnsatzFunctionType;
    using Jacobian               = typename Traits::Jacobian;
    using JacobianColType        = typename Traits::JacobianColType;
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;
    static constexpr int gridDim = Traits::gridDim;

    template <typename WrtType, int I>
    static constexpr bool hasCoeff = DerivativeDirections::HasCoeff<WrtType, gridDim, I>;

    /** \brief Return the function value at the i-th bound integration point*/
    FunctionReturnType evaluateFunction(long unsigned i) {
      const auto& N = impl().basis.getFunction(i);
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
    requires(hasCoeff<Wrt<Args...>, 2>and DerivativeDirections::HasNoSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                             TransformWith<TransformArgs...>&& transArgs,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);

      AnsatzFunctionJacobian dN;
      // Check if a matrix is given to transform derivatives. Otherwise we do nothing
      if constexpr (sizeof...(TransformArgs) > 0)
        dN = dNraw * std::get<0>(transArgs.args);
      else
        dN = dNraw;
      return impl().evaluateSecondDerivativeWRTCoeffs(N, dN, std::get<0>(along.args), coeffsIndices.args);
    }

    template <typename... Args, typename... AlongArgs, typename... Indices, typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 2>and DerivativeDirections::HasNoSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      return evaluateDerivative(localOrIpId, std::forward<Wrt<Args...>>(args), std::forward<Along<AlongArgs...>>(along),
                                transformWith(), std::forward<CoeffIndices<Indices...>>(coeffsIndices));
    }
    //
    /** \brief Function to forward the call of one spatial derivative and two derivative wrt. coefficients.
     * You have to pass a along argument which specifies the direction wher this derivative is applied */
    template <typename... Args, typename... TransformArgs, typename... AlongArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 2>and DerivativeDirections::HasOneSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                             TransformWith<TransformArgs...>&& transArgs,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      transformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));
      constexpr DerivativeDirections::Counter<gridDim> counter
          = DerivativeDirections::countInTuple<Wrt<Args...>, gridDim>();
      constexpr int spatialIndex = DerivativeDirections::findSingleSpatial(counter.spatialDerivatives);
      if constexpr (spatialIndex < gridDim)
        return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
            N, dNTransformed, std::get<0>(along.args), coeffsIndices.args, spatialIndex);
      else if constexpr (counter.coeffDerivatives == 2 and counter.spatialall == 1)
        return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(N, dNTransformed, std::get<0>(along.args),
                                                                             coeffsIndices.args);
    }

    template <typename... TransformArgs>
    void transformDerivatives(const AnsatzFunctionJacobian& dNraw, TransformWith<TransformArgs...>&& transArgs) const {
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
        const AnsatzFunctionJacobian& dN = impl().basis.getJacobian(localOrIpId);
        const AnsatzFunctionType& N      = impl().basis.getFunction(localOrIpId);
        return std::make_tuple(std::ref(N), std::ref(dN));
      } else
        static_assert(std::is_same_v<DomainTypeOrIntegrationPointIndex,
                                     DomainType> or std::is_same_v<DomainTypeOrIntegrationPointIndex, int>,
                      "The argument you passed should be an id for the integration point or the point where the "
                      "derivative should be evaluated");
    }

    /** \brief Function to forward the call of one spatial derivative and no derivative wrt. coefficients  */
    template <typename... Args, typename... TransformArgs, typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 0>and DerivativeDirections::HasOneSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             TransformWith<TransformArgs...>&& transArgs) const {
      const auto& [N, dN]    = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      constexpr auto counter = DerivativeDirections::countInTuple<Wrt<Args...>, gridDim>();
      // Check if a matrix is given to transform derivatives. Otherwise we do nothing
      transformDerivatives(dN, std::forward<TransformWith<TransformArgs...>>(transArgs));
      constexpr auto singleSpatial = DerivativeDirections::findSingleSpatial<gridDim>(counter.spatialDerivatives);
      if constexpr (counter.spatialall == 1)
        return impl().evaluateDerivativeWRTSpaceAllImpl(N, dNTransformed);
      else if constexpr (singleSpatial < gridDim)
        return impl().evaluateDerivativeWRTSpaceSingleImpl(N, dNTransformed, singleSpatial);
      else
        static_assert((singleSpatial < gridDim) or (counter.spatialall == 1),
                      "This currently only supports first order spatial derivatives");
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

    /** \brief Function to forward the call of one spatial derivative and one derivative wrt. coefficients  */
    template <typename... Args, typename... TransformArgs, typename... Indices,
              typename DomainTypeOrIntegrationPointIndex>
    requires(hasCoeff<Wrt<Args...>, 1>and DerivativeDirections::HasOneSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args,
                                                             TransformWith<TransformArgs...>&& transArgs,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      transformDerivatives(dNraw, std::forward<TransformWith<TransformArgs...>>(transArgs));
      constexpr DerivativeDirections::Counter<gridDim> counter
          = DerivativeDirections::countInTuple<Wrt<Args...>, gridDim>();
      constexpr int spatialIndex = DerivativeDirections::findSingleSpatial(counter.spatialDerivatives);
      if constexpr (spatialIndex < gridDim)
        return impl().evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(N, dNTransformed, coeffsIndices.args[0],
                                                                               spatialIndex);
      else if constexpr (counter.coeffDerivatives == 1 and counter.spatialall == 1)
        return impl().evaluateDerivativeWRTCoeffsANDSpatialImpl(N, dNTransformed, coeffsIndices.args[0]);
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
    mutable AnsatzFunctionJacobian dNTransformed;

    LocalFunctionImpl const& impl() const  // CRTP
    {
      return static_cast<LocalFunctionImpl const&>(*this);
    }
  };
}  // namespace Ikarus
