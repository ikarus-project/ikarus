//
// Created by alex on 3/17/22.
//

#pragma once
#include <concepts>
namespace Ikarus {
  namespace DerivativeDirections {
    static struct Coeffs {
    } coeffs;
    static struct Spatial { } spatial; }  // namespace DerivativeDirections

  template <typename... Args>
  struct CoeffIndices {
    std::array<std::common_type_t<std::remove_cvref_t<Args>...>, sizeof...(Args)> args;
  };

  template <typename... Args>
  auto coeffIndices(Args&&... args) {
    return CoeffIndices<Args&&...>({std::forward<Args>(args)...});
  }

  template <typename... Args>
  struct Wrt {
    std::tuple<std::remove_cvref_t<Args>...> args;
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

  struct Counter {
    int coeffDerivatives;
    int spatialDerivatives;
  };

  template <typename WrtType>
  consteval Counter countInTuple() {
    Counter counter{0, 0};
    Dune::Hybrid::forEach(
        Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(WrtType::args)>>()), [&](auto i) {
          using CurrentType = std::tuple_element_t<i, decltype(WrtType::args)>;
          if constexpr (std::is_same_v<DerivativeDirections::Coeffs, CurrentType>)
            ++counter.coeffDerivatives;
          else if constexpr (std::is_same_v<DerivativeDirections::Spatial, CurrentType>)
            ++counter.spatialDerivatives;
          else
            static_assert(std::is_same_v<DerivativeDirections::Coeffs,
                                         CurrentType> or std::is_same_v<DerivativeDirections::Spatial, CurrentType>,
                          "You are only allowed to pass coeffs and spatial struct");
        });
    return counter;
  }

  template <typename WrtType,int I>
  concept HasCoeff =  (countInTuple<WrtType>().coeffDerivatives==I);
  template <typename WrtType,int I>
  concept HasSpatial =  (countInTuple<WrtType>().spatialDerivatives==I);

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
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

    FunctionReturnType evaluateFunction(long unsigned i) {
      const auto& N = impl().basis.getFunction(i);
      return impl().evaluateFunctionImpl(N);
    }

    FunctionReturnType evaluateFunction(const DomainType& local) {
      AnsatzFunctionType N;
      impl().basis.evaluateFunction(local, N);
      return impl().evaluateFunctionImpl(N);
    }

    template <typename... Args, typename... TransformArgs, typename... Indices>
    requires(HasCoeff<Wrt<Args...>,1> and HasSpatial<Wrt<Args...>,1>)
    auto evaluateDerivative(long unsigned gpIndex, const FunctionReturnType& val,
                                                 Wrt<Args...>&& args, TransformWith<TransformArgs...>&& transArgs,
                                                 CoeffIndices<Indices...>&& coeffsIndices) const {
      constexpr Counter counter = countInTuple<Wrt<Args...>>();
      if constexpr (counter.coeffDerivatives == 0 and counter.spatialDerivatives == 1)
        return evaluateDerivative(gpIndex, val, transArgs.args[0]);
      else if constexpr (counter.coeffDerivatives == 1 and counter.spatialDerivatives == 1)
        return impl().evaluateDerivativeWRTCoeffsANDSpatialImpl(gpIndex, val, std::get<0>(transArgs.args),
                                                                coeffsIndices.args[0]);
    }

    template <typename... Args, typename... Indices>
    requires(HasCoeff<Wrt<Args...>,1>
             and HasSpatial<Wrt<Args...>,0>) auto evaluateDerivative(long unsigned gpIndex, const FunctionReturnType& val,
                                                 Wrt<Args...>&& args, CoeffIndices<Indices...>&& coeffsIndices) const {
      return impl().evaluateDerivativeWRTCoeffsImpl(gpIndex, val, coeffsIndices.args[0]);
    }

    template <typename... Args, typename... AlongArgs, typename... Indices>
    requires((countInTuple<Wrt<Args...>>().coeffDerivatives == 2)
             and (countInTuple<Wrt<Args...>>().spatialDerivatives
                  == 0)) auto evaluateDerivative(long unsigned gpIndex, const FunctionReturnType& val,
                                                 Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                 CoeffIndices<Indices...>&& coeffsIndices) const {
      return impl().evaluateSecondDerivativeWRTCoeffs(gpIndex, val, std::get<0>(along.args), coeffsIndices.args);
    }

    template <typename... Args, typename... TransformArgs, typename... AlongArgs, typename... Indices>
    requires((countInTuple<Wrt<Args...>>().coeffDerivatives == 2)
             and (countInTuple<Wrt<Args...>>().spatialDerivatives
                  == 1)) auto evaluateDerivative(long unsigned gpIndex, const FunctionReturnType& val,
                                                 Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                                 TransformWith<TransformArgs...>&& transArgs,
                                                 CoeffIndices<Indices...>&& coeffsIndices) const {
      return impl().evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
          gpIndex, val, std::get<0>(along.args), std::get<0>(transArgs.args), coeffsIndices.args);
    }

    template <typename... Args, typename... TransformArgs>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives == 0) Jacobian
        evaluateDerivative(long unsigned gpIndex, const FunctionReturnType& val, Wrt<Args...>&& args,
                           TransformWith<TransformArgs...>&& transArgs)
    const {
      const AnsatzFunctionJacobian dN = (impl().basis.getJacobian(gpIndex) * std::get<0>(transArgs.args));
      constexpr Counter counter       = countInTuple<Wrt<Args...>>();

      if constexpr (counter.spatialDerivatives == 1)
        return impl().evaluateDerivativeWRTSpaceImpl(impl().basis.getFunction(gpIndex), dN, val);
      else
        static_assert(counter.spatialDerivatives == 1, "This currently only supports first order spatial derivatives");
    }

    template <typename... Args, typename... TransformArgs>
    requires((countInTuple<Wrt<Args...>>().coeffDerivatives == 0)
             and (countInTuple<Wrt<Args...>>().spatialDerivatives == 1)) Jacobian
        evaluateDerivative(long unsigned gpIndex, Wrt<Args...>&& args, TransformWith<TransformArgs...>&& transArgs)
    const {
      const FunctionReturnType val = impl().evaluateFunctionImpl(impl().basis.getFunction(gpIndex));
      return evaluateDerivative(gpIndex, val, std::forward<Wrt<Args...>>(args),
                                std::forward<TransformWith<TransformArgs...>>(transArgs));
    }

    template <typename... Args, typename... Indices>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives
             == 1) auto evaluateDerivative(long unsigned gpIndex, Wrt<Args...>&& args,
                                           CoeffIndices<Indices...>&& coeffsIndices) const {
      const FunctionReturnType val = impl().evaluateFunctionImpl(impl().basis.getFunction(gpIndex));
      return evaluateDerivative(gpIndex, val, std::forward<Wrt<Args...>>(args),
                                std::forward<CoeffIndices<Indices...>>(coeffsIndices));
    }

    template <typename... Args, typename... AlongArgs, typename... Indices>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives
             == 2) auto evaluateDerivative(long unsigned gpIndex, Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                           CoeffIndices<Indices...>&& coeffsIndices) const {
      FunctionReturnType val = impl().evaluateFunctionImpl(impl().basis.getFunction(gpIndex));
      return evaluateDerivative(gpIndex, val, std::forward<Wrt<Args...>>(args),
                                std::forward<Along<AlongArgs...>>(along),
                                std::forward<CoeffIndices<Indices...>>(coeffsIndices));
    }

    auto viewOverIntegrationPoints() { return impl().basis.viewOverIntegrationPoints(); }

  private:
    LocalFunctionImpl const& impl() const  // CRTP
    {
      return static_cast<LocalFunctionImpl const&>(*this);
    }
  };
}  // namespace Ikarus
