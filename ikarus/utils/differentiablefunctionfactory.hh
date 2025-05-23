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
    auto asPtr = toSharedPointer(std::forward<Assembler>(as));

    using FERequirement = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;

    [[maybe_unused]] auto KFunction = [dbcOption, assembler = asPtr, affordances](const FERequirement& p) -> auto& {
      return assembler->matrix(p, affordances.matrixAffordance(), dbcOption);
    };

    [[maybe_unused]] auto residualFunction = [dbcOption, assembler = asPtr,
                                              affordances](const FERequirement& p) -> auto& {
      return assembler->vector(p, affordances.vectorAffordance(), dbcOption);
    };

    auto reqArg = FERequirement(); // This is only created to help  makeDifferentiableFunction deduce the argument type
                                   // of the functions. therefore, no valid values needed

    if constexpr (affordances.hasScalarAffordance) {
      [[maybe_unused]] auto energyFunction = [assembler = asPtr, affordances](const FERequirement& p) -> auto& {
        return assembler->scalar(p, affordances.scalarAffordance());
      };

      return makeDifferentiableFunction(functions(energyFunction, residualFunction, KFunction), reqArg);
    } else
      return makeDifferentiableFunction(functions(residualFunction, KFunction), reqArg);
  }

  template <typename Assembler>
  static auto op(Assembler&& as, DBCOption dbcOption) {
    auto asPtr = toSharedPointer(std::forward<Assembler>(as));
    if (as->boundToAffordanceCollection())
      return op(asPtr, as->affordanceCollection(), dbcOption);

    DUNE_THROW(Dune::InvalidStateException,
               "Assembler has to be bound to an affordance collection before you can call "
               "this method");
  }

  template <typename Assembler>
  static auto op(Assembler&& as) {
    auto asPtr = toSharedPointer(std::forward<Assembler>(as));

    if (asPtr->boundToAffordanceCollection() and asPtr->boundToDBCOption())
      return op(asPtr, asPtr->affordanceCollection(), asPtr->dBCOption());
    DUNE_THROW(Dune::InvalidStateException,
               "Assembler has to be bound to an affordance collection and to an "
               "DBCOption before you can call "
               "this method");
  }

private:
  // Convert to shared pointer if needed
  template <typename Assembler>
  static decltype(auto) toSharedPointer(Assembler&& as) {
    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> ||
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      return std::forward<Assembler>(as);
    } else {
      return std::make_shared<std::remove_cvref_t<Assembler>>(std::forward<Assembler>(as));
    }
  }
};

} // namespace Ikarus
