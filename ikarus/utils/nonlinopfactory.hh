// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinopfactory.hh
 * \brief Contains the generic NonLinearOperatorFactory class.
 */

#pragma once

#include <utility>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/nonlinearoperator.hh>

namespace Ikarus {

struct NonLinearOperatorFactory
{
  template <typename Assembler, typename... Affordances>
  static auto op(Assembler&& as, typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement& req,
                 AffordanceCollection<Affordances...> affordances, DBCOption dbcOption) {
    auto assemblerPtr = [as]() {
      if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                    traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value)
        return as;
      else
        return std::make_shared<std::remove_cvref_t<Assembler>>(std::forward<Assembler>(as));
    }();

    using FERequirement             = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;
    [[maybe_unused]] auto KFunction = [dbcOption, assembler = assemblerPtr, affordances](
                                          typename FERequirement::SolutionVectorType& globalSol,
                                          typename FERequirement::ParameterType& parameter) -> auto& {
      FERequirement req;
      req.insertGlobalSolution(globalSol).insertParameter(parameter);

      return assembler->matrix(req, affordances.matrixAffordance(), dbcOption);
    };

    [[maybe_unused]] auto residualFunction = [dbcOption, assembler = assemblerPtr, affordances](
                                                 typename FERequirement::SolutionVectorType& globalSol,
                                                 typename FERequirement::ParameterType& parameter) -> auto& {
      FERequirement req;
      req.insertGlobalSolution(globalSol).insertParameter(parameter);
      return assembler->vector(req, affordances.vectorAffordance(), dbcOption);
    };

    assert(req.populated() && " Before you calls this method you have to pass populated fe requirements");
    if constexpr (affordances.hasScalarAffordance) {
      [[maybe_unused]] auto energyFunction = [assembler = assemblerPtr, affordances](
                                                 typename FERequirement::SolutionVectorType& globalSol,
                                                 typename FERequirement::ParameterType& parameter) -> auto& {
        FERequirement req;
        req.insertGlobalSolution(globalSol).insertParameter(parameter);

        return assembler->scalar(req, affordances.scalarAffordance());
      };
      return NonLinearOperator(functions(std::move(energyFunction), std::move(residualFunction), std::move(KFunction)),
                               parameter(req.globalSolution(), req.parameter()));
    } else
      return NonLinearOperator(functions(std::move(residualFunction), std::move(KFunction)),
                               parameter(req.globalSolution(), req.parameter()));
  }

  template <typename Assembler>
  static auto op(Assembler&& as, DBCOption dbcOption) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to a fe requirement and an affordance collection before you can call "
                 "this method");
    };
    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (as->boundToRequirement() and as->boundToAffordanceCollection()) {
        return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), dbcOption);
      } else {
        ex();
      }
    } else {
      if (as->boundToRequirement() and as->boundToAffordanceCollection()) {
        return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), dbcOption);
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
                 "Assembler has to be bound to a fe requirement to an affordance collection and to an "
                 "DBCOption before you can call "
                 "this method");
    };
    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (not as->bound())
        ex();
      return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), as->dBCOption());
    } else {
      if (not as.bound())
        ex();
      return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), as->dBCOption());
    }
  }

  template <typename Assembler, typename... Affordances>
  static auto op(Assembler&& as, AffordanceCollection<Affordances...> affordances,
                 DBCOption dbcOption = DBCOption::Full) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to a fe requirement before you can call "
                 "this method");
    };

    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (not as->boundToRequirement())
        ex();
      return op(std::forward<Assembler>(as), as->requirement(), affordances, dbcOption);
    } else {
      if (not as.boundToRequirement())
        ex();
      return op(std::forward<Assembler>(as), as.requirement(), affordances, dbcOption);
    }
  }

  template <typename Assembler>
  static auto op(Assembler&& as, typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement& req,
                 DBCOption dbcOption) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to an affordance collection before you can call "
                 "this method");
    };

    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (not as->boundToAffordanceCollection())
        ex();
      return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), dbcOption);
    } else {
      if (not as.boundToAffordanceCollection())
        ex();
      return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), dbcOption);
    }
  }
};

} // namespace Ikarus
