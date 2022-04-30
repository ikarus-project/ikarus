//
// Created by Alex on 26.04.2022.
//

#pragma once
#include <cstddef>

#include <ikarus/utils/traits.hh>
#include <dune/typetree/typetree.hh>
#include <dune/istl/multitypeblockvector.hh>
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

  struct DerivativeNoOp{};

    [[maybe_unused]] static struct SpatialAll { } spatialAll;

    struct SpatialPartial {
      size_t index{};
    };


  template<std::size_t I>
  struct SingleCoeff {
    Dune::MultiTypeBlockVector<decltype(Dune::TypeTree::treePath(std::declval<Dune::index_constant<I>>(),std::declval<size_t>()))> index{};
  };



  template<std::size_t I,std::size_t J>
  struct TwoCoeff {
    Dune::MultiTypeBlockVector<decltype(Dune::TypeTree::treePath(std::declval<Dune::index_constant<I>>(),std::declval<size_t>())),
               decltype(Dune::TypeTree::treePath(std::declval<Dune::index_constant<J>>(),std::declval<size_t>()))> index{};
  };

    SpatialPartial spatial(size_t i);

  template<std::size_t I>
  SingleCoeff<I> coeff(Dune::index_constant<I> iObj, size_t i)
  {
    using namespace Dune::Indices;
    SingleCoeff<I> coeffs;
    std::get<1>(coeffs.index[_0]._data) = i;
    return coeffs;
  }
  template<std::size_t I,std::size_t J>
  TwoCoeff<I,J> coeff(Dune::index_constant<I> iObj, size_t i,Dune::index_constant<J> jObj, size_t j)
  {
    using namespace Dune::Indices;
    TwoCoeff<I,J> coeffs;
    std::get<1>(coeffs.index[_0]._data) = i;
    std::get<1>(coeffs.index[_1]._data) = j;
    return coeffs;
  }

  SingleCoeff<0> coeff( size_t i);
  TwoCoeff<0,0> coeff( size_t i, size_t j);

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

    template<typename Type>
    concept isSpatial =
       std::is_same_v<Type,DerivativeDirections::SpatialPartial> or
          std::is_same_v<Type,DerivativeDirections::SpatialAll>;

  template<typename Type>
  concept isCoeff = Std::IsSpecializationNonTypes<SingleCoeff, Type>::value or
      Std::IsSpecializationNonTypes<TwoCoeff, Type>::value;


    struct ConstExprCounter {
      int singleCoeffDerivs{};
      int twoCoeffDerivs{};
      int spatialDerivs{};
      int spatialAllCounter{};

      consteval int orderOfDerivative() const {
        return singleCoeffDerivs + 2 * twoCoeffDerivs + spatialDerivs + spatialAllCounter;
      }
    };

    template <typename WrtType>
    consteval ConstExprCounter countDerivativesType() {
      ConstExprCounter counter{};
      using Tuple               = typename WrtType::Args;
      counter.singleCoeffDerivs = Ikarus::Std::countTypeSpecialization<SingleCoeff,Tuple>();
      counter.twoCoeffDerivs    = Ikarus::Std::countTypeSpecialization<TwoCoeff,Tuple>();
      counter.spatialDerivs     = Ikarus::Std::countType<Tuple, SpatialPartial>();
      counter.spatialAllCounter        = Ikarus::Std::countType<Tuple, SpatialAll>();
      return counter;
    }

    template <typename WrtType>
    auto extractCoeffIndices(WrtType&& wrt) {
      if constexpr (Std::hasTypeSpecialization<SingleCoeff, typename std::remove_reference_t<WrtType>::Args>())
        return Std::getSpecialization<SingleCoeff>(wrt.args).index;  // returns single int
      else if constexpr (Std::hasTypeSpecialization<TwoCoeff, typename std::remove_reference_t<WrtType>::Args>())
        return Std::getSpecialization<TwoCoeff>(wrt.args).index;  // return std::array<size_t,2>
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
        = (countDerivativesType<WrtType>().spatialDerivs == 0 and countDerivativesType<WrtType>().spatialAllCounter == 0);

    template <typename WrtType>
    concept HasOneSpatialAll = countDerivativesType<WrtType>()
    .spatialAllCounter == 1;

    template <typename WrtType>
    concept HasOneSpatialSingle = countDerivativesType<WrtType>()
    .spatialDerivs == 1;

    template <typename WrtType>
    concept HasOneSpatial = HasOneSpatialSingle<WrtType> or HasOneSpatialAll<WrtType>;

  }  // namespace DerivativeDirections

    template <typename LF>
  concept IsUnaryExpr = LF::children==1;

  template <typename LF>
  concept IsBinaryExpr = LF::children==2;

  using Arithmetic = Dune::index_constant<100>;
  static constexpr auto arithmetic = Arithmetic{};

template <typename LF>
concept IsArithmeticExpr = LF::Ids::value==arithmetic;

template <typename E1,typename E2>
class LocalFunctionScale;

template <typename LocalFunctionImpl>
class LocalFunctionInterface;

template <typename LocalFunctionImpl>
concept LocalFunction = requires
{
  typename LocalFunctionImpl::Traits;
  LocalFunctionImpl::Traits::valueSize;
  typename LocalFunctionImpl::Traits::DomainType;
  typename LocalFunctionImpl::Ids;
};

template <typename LF>
concept IsScaleExpr = Std::isSpecialization<LocalFunctionScale,LF>::value;

template <typename LF>
concept IsNonArithmeticLeafNode = LF::isLeaf==true and !IsArithmeticExpr<LF>;

template <typename... LF>
concept IsLocalFunction = (LocalFunction<LF> and ...);


}  // namespace Ikarus