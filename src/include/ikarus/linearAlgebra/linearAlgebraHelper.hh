//
// Created by alex on 3/14/22.
//

#pragma once

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <Eigen/Core>

#include <ikarus/utils/concepts.hh>

namespace Ikarus {

  /* Returns the total correction size of a block vector with a Manifold as type */
  template <typename Type>
  size_t correctionSize(const Dune::BlockVector<Type>& a) requires requires {
    Type::correctionSize;
  }
  { return a.size() * Type::correctionSize; }

  /* Enables the += operator for Dune::BlockVector += Eigen::Vector */
  template <typename Type, typename Derived>
  Dune::BlockVector<Type>& operator+=(Dune::BlockVector<Type>& a, const Eigen::MatrixBase<Derived>& b) requires(
      Ikarus::Concepts::AddAssignAble<Type, decltype(b.template segment<Type::correctionSize>(0))>and requires() {
        Type::correctionSize;
      }) {
    for (auto i = 0U; i < a.size(); ++i)
      a[i] += b.template segment<Type::correctionSize>(i * Type::correctionSize);
    return a;
  }

  /* Enables the += operator for une::MultiTypeBlockVector += Eigen::Vector */
  template <typename... Types, typename Derived>
  Dune::MultiTypeBlockVector<Types...>& operator+=(Dune::MultiTypeBlockVector<Types...>& a,
                                                   const Eigen::MatrixBase<Derived>& b) {
    using namespace Dune::Indices;
    size_t posStart = 0;
    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<a.size()>()), [&](const auto i) {
      const size_t size = correctionSize(a[i]);
      a[i] += b(Eigen::seqN(posStart, size));
      posStart += size;
    });

    return a;
  }

}  // namespace Ikarus