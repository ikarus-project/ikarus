//
// Created by alex on 3/17/22.
//

#pragma once
#include <concepts>

#include "ikarus/LocalBasis/localBasis.h"
namespace Ikarus {
  namespace DerivativeDirections {
    constexpr std::integral_constant<int, -1> coeffs = {};

    //    template <int i>
    //    requires(i >= 0 and i < 3) using spatial             = std::integral_constant<int, i>;
    //    constexpr spatial<0> spatial0                        = {};
    //    constexpr spatial<1> spatial1                        = {};
    //    constexpr spatial<2> spatial2                        = {};
    constexpr std::integral_constant<int, -2> spatialall = {};

    struct spatialPartial {
      using value = std::integral_constant<int, -3>;
      int index{};
    };

    spatialPartial spatial(int i) { return {i}; }

    template <int Dim>
    struct Counter {
      int coeffDerivatives{};
      std::array<int, Dim> spatialDerivatives{};
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
    Counter<gridDim> countDynamicSpatialDerivativesInTuple(const WrtType& wrt) {
      Counter<gridDim> counter{};
      Dune::Hybrid::forEach(
          Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(WrtType::args)>>()), [&](auto i) {
            using currentDerivType = std::tuple_element_t<i, typename WrtType::Args>;
            if constexpr (std::is_same_v<currentDerivType, spatialPartial>)
              ++counter.spatialDerivatives[std::get<i>(wrt.args).index];
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
    concept HasOneSpatial
        = (countInTuple<WrtType, gridDim>().dynamicspatial == 1 or countInTuple<WrtType, gridDim>().spatialall == 1);

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
    using ctype                  = typename Traits::ctype;
    using FunctionReturnType     = typename Traits::FunctionReturnType;
    using AnsatzFunctionType     = typename Traits::AnsatzFunctionType;
    using Jacobian               = typename Traits::Jacobian;
    using JacobianColType        = typename Traits::JacobianColType;
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;
    using CoeffDerivMatrix       = typename Traits::CoeffDerivMatrix;
    /** \brief Matrix to transform the ansatz function Jacobian to world coordinates*/

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
    requires(hasCoeff<Wrt<Args...>, 2>and DerivativeDirections::HasNoSpatial<
             Wrt<Args...>, gridDim>) auto evaluateDerivative(const DomainTypeOrIntegrationPointIndex& localOrIpId,
                                                             Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                             TransformWith<TransformArgs...>&& transArgs,
                                                             CoeffIndices<Indices...>&& coeffsIndices) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);

      // Check if a matrix is given to transform derivatives. Otherwise we do nothing
      if constexpr (sizeof...(TransformArgs) > 0) {
        AnsatzFunctionJacobian dN = (dNraw * std::get<0>(transArgs.args)).eval();
        return tryCallSecondDerivativeWRTCoeffs(N, dN, std::get<0>(along.args), coeffsIndices.args);
      } else
        return tryCallSecondDerivativeWRTCoeffs(N, dNraw, std::get<0>(along.args), coeffsIndices.args);
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
      const DerivativeDirections::Counter<gridDim> counter
          = DerivativeDirections::countDynamicSpatialDerivativesInTuple<Wrt<Args...>, gridDim>(
              std::forward<Wrt<Args...>>(args));
      constexpr DerivativeDirections::Counter<gridDim> counterExpr
          = DerivativeDirections::countInTuple<Wrt<Args...>, gridDim>();
      const int spatialIndex = DerivativeDirections::findSingleSpatial(counter.spatialDerivatives);
      if (spatialIndex < gridDim)
        return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
            N, dNTransformed, std::get<0>(along.args), coeffsIndices.args, spatialIndex);
      else if constexpr (counterExpr.coeffDerivatives == 2 and counterExpr.spatialall == 1)
        return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(N, dNTransformed, std::get<0>(along.args),
                                                                             coeffsIndices.args);
      else
        __builtin_unreachable();
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
        const AnsatzFunctionJacobian& dN = impl().basis.evaluateJacobian(localOrIpId);
        const AnsatzFunctionType& N      = impl().basis.evaluateFunction(localOrIpId);
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
      const auto& [N, dN] = evaluateFunctionAndDerivativeWithIPorCoord(localOrIpId);
      const DerivativeDirections::Counter<gridDim> counter
          = DerivativeDirections::countDynamicSpatialDerivativesInTuple<Wrt<Args...>, gridDim>(
              std::forward<Wrt<Args...>>(args));
      constexpr DerivativeDirections::Counter<gridDim> counterConstExpr
          = DerivativeDirections::countInTuple<Wrt<Args...>, gridDim>();
      // Check if a matrix is given to transform derivatives. Otherwise we do nothing
      transformDerivatives(dN, std::forward<TransformWith<TransformArgs...>>(transArgs));
      const auto singleSpatial = DerivativeDirections::findSingleSpatial<gridDim>(counter.spatialDerivatives);
      if constexpr (counterConstExpr.spatialall == 1)
        return impl().evaluateDerivativeWRTSpaceAllImpl(N, dNTransformed);
      else if (singleSpatial < gridDim)
        return impl().evaluateDerivativeWRTSpaceSingleImpl(N, dNTransformed, singleSpatial);
      else
        static_assert((singleSpatial < gridDim) or (counterConstExpr.spatialall == 1),
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
      const DerivativeDirections::Counter<gridDim> counter
          = DerivativeDirections::countDynamicSpatialDerivativesInTuple<Wrt<Args...>, gridDim>(
              std::forward<Wrt<Args...>>(args));
      constexpr DerivativeDirections::Counter<gridDim> counterExpr
          = DerivativeDirections::countInTuple<Wrt<Args...>, gridDim>();
      const int spatialIndex = DerivativeDirections::findSingleSpatial(counter.spatialDerivatives);
      if constexpr (counterExpr.coeffDerivatives == 1 and counterExpr.spatialall != 1)
        return impl().evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(N, dNTransformed, coeffsIndices.args[0],
                                                                      spatialIndex);
      else if constexpr (counterExpr.coeffDerivatives == 1 and counterExpr.spatialall == 1)
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
