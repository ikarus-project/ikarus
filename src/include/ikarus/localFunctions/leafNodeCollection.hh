// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include "meta.hh"

#include <dune/common/indices.hh>

#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/utils/traits.hh>
namespace Ikarus {
  namespace Impl {

    template <typename LocalFunctionImpl>
    class LocalFunctionInterface;

    template <typename LF>
    requires(!std::is_arithmetic_v<LF>) consteval int countNonArithmeticLeafNodesImpl() {
      if constexpr (Std::isSpecialization<std::tuple, typename LF::Ids>::value) {
        constexpr auto predicate = []<typename Type>(Type) { return Type::value != Ikarus::arithmetic; };
        return std::tuple_size_v<decltype(Std::filter(typename LF::Ids(), predicate))>;
      } else
        return 1;
    }

    template <typename LF>
    requires LocalFunction<LF>
    auto collectNonArithmeticLeafNodesImpl(LF&& a) {
      using LFRaw = std::remove_cvref_t<LF>;
      //    static_assert(LocalFunction<LF>,"Only passing LocalFunctions allowed");
      if constexpr (IsBinaryExpr<LFRaw>)
        return std::tuple_cat(collectNonArithmeticLeafNodesImpl(a.l()), collectNonArithmeticLeafNodesImpl(a.r()));
      else if constexpr (IsUnaryExpr<LFRaw>)
        return std::make_tuple(collectNonArithmeticLeafNodesImpl(a.m()));
      else if constexpr (IsArithmeticExpr<LFRaw>)
        return std::make_tuple();
      else if constexpr (IsNonArithmeticLeafNode<LFRaw>)
        return std::make_tuple(std::ref(a));
      else
        static_assert("There are currently no other expressions. Thus you should not end up here.");
    }

  }  // namespace Impl

  template <typename LF>
  consteval int countNonArithmeticLeafNodes(const LocalFunctionInterface<LF>&) {
    return Impl::countNonArithmeticLeafNodesImpl<LF>();
  }

  template <typename LF>
  requires LocalFunction<LF>
  auto collectNonArithmeticLeafNodes(LF&& a) {
    return Std::makeNestedTupleFlatAndStoreReferences(Impl::collectNonArithmeticLeafNodesImpl(a.impl()));
  }

  /** This class contains the collection of leaf nodes of a local function expression */
  template <typename LF>
  requires LocalFunction<LF>
  struct LocalFunctionLeafNodeCollection {
    using LFRaw = std::remove_cvref_t<LF>;
    /* Since we need to enable perfect forwarding we have to implement this universal constructor. We also constrain it
     * with requires LocalFunction<LF_> to only allow it for local function types. Without this template and a signature
     * as LocalFunctionLeafNodeCollection( LF&& lf): ... perfect forwarding is not working for constructors. See
     * https://eel.is/c++draft/temp.deduct.call#3
     * */
    template <typename LF_>
    requires LocalFunction<LF_> LocalFunctionLeafNodeCollection(LF_&& lf)
        : leafNodes{collectNonArithmeticLeafNodes(std::forward<LF_>(lf))} {}

    template <std::size_t I = 0>
    requires(Std::countType<typename LFRaw::Ids, Dune::index_constant<I>>()
             == 1) auto& coefficientsRef(Dune::index_constant<I> = Dune::index_constant<I>()) {
      static_assert(
          Std::countType<typename LFRaw::Ids, Dune::index_constant<I>>() == 1,
          "Non-const coefficientsRef() can only be called, if there is only one node with the given leaf node ID.");
      return std::get<I>(leafNodes).coefficientsRef();
    }

    template <std::size_t I = 0>
    const auto& coefficientsRef(Dune::index_constant<I> = Dune::index_constant<0UL>()) const {
      return std::get<I>(leafNodes).coefficientsRef();
    }

    template <typename Derived, std::size_t I = 0>
    void addToCoeffs(const Eigen::MatrixBase<Derived>& correction,
                     Dune::index_constant<I> = Dune::index_constant<0UL>()) {
      Dune::Hybrid::forEach(leafNodes, [&]<typename LFI>(LFI& lfi) {
        if constexpr (LFI::Ids::value == I) lfi.coefficientsRef() += correction;
      });
    }

    // This function updates the embedding values of the underlying coefficients
    template <typename Derived, std::size_t I = 0>
    void addToCoeffsInEmbedding(const Eigen::MatrixBase<Derived>& correction,
                                Dune::index_constant<I> = Dune::index_constant<0UL>()) {
      Dune::Hybrid::forEach(leafNodes, [&]<typename LFI>(LFI& lfi) {
        if constexpr (LFI::Ids::value == I) addInEmbedding(lfi.coefficientsRef(), correction);
      });
    }

    template <std::size_t I = 0>
    auto& basis(Dune::index_constant<I> = Dune::index_constant<0UL>()) {
      return std::get<I>(leafNodes).basis();
    }

    template <std::size_t I = 0>
    auto& node(Dune::index_constant<I> = Dune::index_constant<0UL>()) {
      return std::get<I>(leafNodes);
    }

    constexpr auto size() { return Dune::index_constant<std::tuple_size_v<LeafNodeTuple>>(); }

    using LeafNodeTuple = decltype(collectNonArithmeticLeafNodes(std::declval<LF&&>()));

  private:
    LeafNodeTuple leafNodes;
  };

  template <typename LF>
  requires LocalFunction<LF>
  auto collectLeafNodeLocalFunctions(LF&& lf) { return LocalFunctionLeafNodeCollection<LF>(std::forward<LF>(lf)); }
}  // namespace Ikarus
