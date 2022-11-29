// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <ikarus/localFunctions/expressions/unaryExpr.hh>
#include <ikarus/localFunctions/meta.hh>
namespace Ikarus {

  template <typename Type>
  requires std::is_arithmetic_v<Type>
  class ConstantExpr : public LocalFunctionInterface<ConstantExpr<Type>> {
  public:
    explicit ConstantExpr(Type val_) : val{val_} {}

    const Type& value() const { return val; }
    Type& value() { return val; }

    auto clone() const { return ConstantExpr(val); }

    template <typename OtherType, size_t ID = 0>
    auto rebindClone(OtherType&&, [[maybe_unused]] Dune::index_constant<ID>&& id = Dune::index_constant<0>()) const {
      if constexpr (Arithmetic::value == ID)
        return ConstantExpr(static_cast<OtherType>(val));
      else
        return clone();
    }

    template <typename OtherType>
    struct Rebind {
      using other = ConstantExpr<OtherType>;
    };

    static constexpr bool isLeaf = true;
    using Ids                    = Arithmetic;
    template <size_t ID_ = 0>
    static constexpr int orderID = ID_ == Arithmetic::value ? linear : constant;

  private:
    Type val;
  };

  template <typename Type>
  struct LocalFunctionTraits<ConstantExpr<Type>> {
    static constexpr int valueSize = 1;
    /** \brief Type for the points for evaluation, usually the integration points */
    using DomainType = Dune::FieldVector<double, 0>;
  };

}  // namespace Ikarus
