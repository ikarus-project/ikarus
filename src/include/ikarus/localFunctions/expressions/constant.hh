/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

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
    auto rebindClone(OtherType&& t, Dune::index_constant<ID>&& id = Dune::index_constant<0>()) const {
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