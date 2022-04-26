//
// Created by Alex on 26.04.2022.
//

#pragma once
#include <cstddef>

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
}  // namespace Ikarus