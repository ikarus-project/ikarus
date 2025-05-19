// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file differentiablefunctionfactory.hh
 * \brief Contains the generic DifferentiableFunctionFactory class.
 */

#pragma once

#include <utility>

#include <dune/functions/common/differentiablefunctionfromcallables.hh>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/derivativetraits.hh>
#include <ikarus/utils/differentiablefunction.hh>

namespace Ikarus {

struct DifferentiableFunctionFactory
{
  template <typename Assembler, typename... Affordances>
  static auto op(Assembler&& as, AffordanceCollection<Affordances...> affordances,
                 DBCOption dbcOption = DBCOption::Full) {
    auto assemblerPtr = [as]() {
      if constexpr (traits::Pointer<std::remove_cvref_t<Assembler>> or
                    traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value)
        return as;
      else
        return std::make_shared<std::remove_cvref_t<Assembler>>(std::forward<Assembler>(as));
    }();

    using FERequirement = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;

    [[maybe_unused]] auto KFunction = [dbcOption, assembler = assemblerPtr,
                                       affordances](const FERequirement& p) -> auto& {
      return assembler->matrix(p, affordances.matrixAffordance(), dbcOption);
    };

    [[maybe_unused]] auto residualFunction = [dbcOption, assembler = assemblerPtr,
                                              affordances](const FERequirement& p) -> auto& {
      return assembler->vector(p, affordances.vectorAffordance(), dbcOption);
    };

    auto reqArg = FERequirement(); // This is only created to help  makeDifferentiableFunction deduce the argument type
                                   // of the functions. therefore, no valid values needed

    if constexpr (affordances.hasScalarAffordance) {
      [[maybe_unused]] auto energyFunction = [assembler = assemblerPtr, affordances](const FERequirement& p) -> auto& {
        return assembler->scalar(p, affordances.scalarAffordance());
      };

      return makeDifferentiableFunction(functions(energyFunction, residualFunction, KFunction), reqArg);
    } else
      return makeDifferentiableFunction(functions(residualFunction, KFunction), reqArg);
  }

  template <typename Assembler>
  static auto op(Assembler&& as, DBCOption dbcOption) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to an affordance collection before you can call "
                 "this method");
    };
    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (as->boundToAffordanceCollection()) {
        return op(std::forward<Assembler>(as), as->affordanceCollection(), dbcOption);
      } else {
        ex();
      }
    } else {
      if (as->boundToAffordanceCollection()) {
        return op(std::forward<Assembler>(as), as.affordanceCollection(), dbcOption);
      } else {
        ex();
      }
    }
    __builtin_unreachable();
  }

  template <typename Assembler>
  static auto op(Assembler&& as) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to an affordance collection and to an "
                 "DBCOption before you can call "
                 "this method");
    };
    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (not(as->boundToAffordanceCollection() and as->boundToDBCOption()))
        ex();
      return op(std::forward<Assembler>(as), as->affordanceCollection(), as->dBCOption());
    } else {
      if (not(as->boundToAffordanceCollection() and as->boundToDBCOption()))
        ex();
      return op(std::forward<Assembler>(as), as.affordanceCollection(), as.dBCOption());
    }
  }
};

} // namespace Ikarus
