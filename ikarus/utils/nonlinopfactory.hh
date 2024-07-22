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
#include <dune/common/integersequence.hh>

namespace Ikarus {

namespace Impl {
  template <typename Check, typename F, typename T>
  T optionalDecider(Check&& c, F&& valueFunction, const std::optional<T>& b) {
    if (c()) {
      return valueFunction();
    } else if (b.has_value()) {
      return b.value();
    } else
      DUNE_THROW(Dune::InvalidStateException,
                 "Neither the assembler is bound to the given type nor was a value of type " + Dune::className<T>() +
                     " passed");
  }
} // namespace Impl

struct NonLinearOperatorFactory
{

  template<typename Assembler, typename... Affordances>
  static auto scalarFunction(std::shared_ptr<Assembler> assemblerPtr,AffordanceCollection<Affordances...> affordances) {
        using FERequirement             = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;

    return [assembler = assemblerPtr, affordances](
                                                 typename FERequirement::SolutionVectorType& globalSol,
                                                 typename FERequirement::ParameterType& parameter) -> auto& {
        FERequirement req;
        req.insertGlobalSolution(globalSol).insertParameter(parameter);

        return assembler->scalar(req, affordances.scalarAffordance());
      };
  }

    template<typename Assembler, typename... Affordances>
  static auto vectorFunction(std::shared_ptr<Assembler> assemblerPtr,AffordanceCollection<Affordances...> affordances,DBCOption dbcOption) {
        using FERequirement             = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;

    return [dbcOption, assembler = assemblerPtr, affordances](
                                                 typename FERequirement::SolutionVectorType& globalSol,
                                                 typename FERequirement::ParameterType& parameter) -> auto& {
      FERequirement req;
      req.insertGlobalSolution(globalSol).insertParameter(parameter);
      return assembler->vector(req, affordances.vectorAffordance(), dbcOption);
    };
  }


    template<typename Assembler, typename... Affordances>
  static auto matrixFunction(std::shared_ptr<Assembler> assemblerPtr,AffordanceCollection<Affordances...> affordances,DBCOption dbcOption) {
        using FERequirement             = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;

    return [dbcOption, assembler = assemblerPtr, affordances](
                                          typename FERequirement::SolutionVectorType& globalSol,
                                          typename FERequirement::ParameterType& parameter) -> auto& {
      FERequirement req;
      req.insertGlobalSolution(globalSol).insertParameter(parameter);

      return assembler->matrix(req, affordances.matrixAffordance(), dbcOption);
    };
  }
//TODO Alex makes this into a singlew function templated with i






  template <size_t... funcs, typename Assembler, typename... Affordances>
  static auto op(Assembler&& as, typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement& req,
                 AffordanceCollection<Affordances...> affordances, DBCOption dbcOption,
                 std::integer_sequence<size_t,funcs...> funcIndices = {}) {
    constexpr int funcIndexSize = funcIndices.size();
    static_assert(Dune::equal(Dune::sorted(funcIndices),funcIndices), "The function indices you request have to be sorted");
    static_assert(funcIndexSize < 4, "The function indices you request have to be less than 4");
    static_assert(Dune::filter([](auto i) { return i < 3; },funcIndices).size() == funcIndexSize,
                  "The function indices you request have to be less than 3");
    using namespace Dune::Indices;

    constexpr bool provideScalar =
        (affordances.hasScalarAffordance and (funcIndexSize == 0 or Dune::contains(funcIndices, _0)));
    constexpr bool provideVector =
        (affordances.hasVectorAffordance and (funcIndexSize == 0 or Dune::contains(funcIndices, _1)));
    constexpr bool provideMatrix =
        (affordances.hasMatrixAffordance and (funcIndexSize == 0 or Dune::contains(funcIndices, _2)));
    auto assemblerPtr = [as]() {
      if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                    traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value)
        return as;
      else
        return std::make_shared<std::remove_cvref_t<Assembler>>(std::forward<Assembler>(as));
    }();

    using FERequirement             = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;
    // [[maybe_unused]] auto KFunction = [dbcOption, assembler = assemblerPtr, affordances](
    //                                       typename FERequirement::SolutionVectorType& globalSol,
    //                                       typename FERequirement::ParameterType& parameter) -> auto& {
    //   FERequirement req;
    //   req.insertGlobalSolution(globalSol).insertParameter(parameter);

    //   return assembler->matrix(req, affordances.matrixAffordance(), dbcOption);
    // };

    // [[maybe_unused]] auto residualFunction = [dbcOption, assembler = assemblerPtr, affordances](
    //                                              typename FERequirement::SolutionVectorType& globalSol,
    //                                              typename FERequirement::ParameterType& parameter) -> auto& {
    //   FERequirement req;
    //   req.insertGlobalSolution(globalSol).insertParameter(parameter);
    //   return assembler->vector(req, affordances.vectorAffordance(), dbcOption);
    // };

    assert(req.populated() && " Before you calls this method you have to pass populated fe requirements");
    if constexpr (provideScalar) {
      // [[maybe_unused]] auto energyFunction = [assembler = assemblerPtr, affordances](
      //                                            typename FERequirement::SolutionVectorType& globalSol,
      //                                            typename FERequirement::ParameterType& parameter) -> auto& {
      //   FERequirement req;
      //   req.insertGlobalSolution(globalSol).insertParameter(parameter);

      //   return assembler->scalar(req, affordances.scalarAffordance());
      // };
      if constexpr (provideVector and provideMatrix)
        return NonLinearOperator(
            functions(scalarFunction(assemblerPtr,affordances), vectorFunction(assemblerPtr,affordances,dbcOption), matrixFunction(assemblerPtr,affordances,dbcOption)),
            parameter(req.globalSolution(), req.parameter()));
      else if constexpr (provideVector)
        return NonLinearOperator(functions(scalarFunction(assemblerPtr,affordances), vectorFunction(assemblerPtr,affordances,dbcOption)),
                                 parameter(req.globalSolution(), req.parameter()));
      else if constexpr (provideMatrix)
        return NonLinearOperator(functions(scalarFunction(assemblerPtr,affordances), matrixFunction(assemblerPtr,affordances,dbcOption)),
                                 parameter(req.globalSolution(), req.parameter()));
      else
        return NonLinearOperator(functions(scalarFunction(assemblerPtr,affordances)),
                                 parameter(req.globalSolution(), req.parameter()));

    } else if constexpr (provideVector) {
      if constexpr (provideMatrix)
        return NonLinearOperator(functions(vectorFunction(assemblerPtr,affordances,dbcOption), matrixFunction(assemblerPtr,affordances,dbcOption)),
                                 parameter(req.globalSolution(), req.parameter()));
      else
        return NonLinearOperator(functions(vectorFunction(assemblerPtr,affordances,dbcOption)),
                                 parameter(req.globalSolution(), req.parameter()));
    } else if constexpr (provideMatrix)
      return NonLinearOperator(functions(matrixFunction(assemblerPtr,affordances,dbcOption)), parameter(req.globalSolution(), req.parameter()));
    else
      static_assert(Dune::AlwaysFalse<Assembler>::value, "You provided unsutiable function indices");
  }

  template <size_t... funcs, typename Assembler, typename... Affordances>
  static auto op(
      Assembler&& as,
      std::optional<
          std::reference_wrapper<typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement>>
          reqArg                                                         = std::nullopt,
      std::optional<AffordanceCollection<Affordances...>> affordancesArg = std::nullopt,
      std::optional<DBCOption> dbCOptionArg                              = std::nullopt,std::integer_sequence<size_t,funcs...> funcIndices = {}) {
    using FERequirement = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;
    auto assemblerPtr   = [as]() {
      if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                    traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value)
        return as;
      else
        return std::make_shared<std::remove_cvref_t<Assembler>>(std::forward<Assembler>(as));
    }();
    FERequirement& req   = Impl::optionalDecider([&]() { return assemblerPtr->boundToRequirement(); },
                                               [&]() mutable { return std::ref(assemblerPtr->requirement()); }, reqArg);
    const auto dbcOption = Impl::optionalDecider([&]() { return assemblerPtr->boundToDBCOption(); },
                                                 [&]() { return assemblerPtr->dBCOption(); }, dbCOptionArg);

    const auto affordances = [&]() {
      using AffoCollectionOfAssembler = typename std::decay_t<decltype(assemblerPtr->affordanceCollection())>;
      if constexpr (sizeof...(Affordances) == 0)
        return assemblerPtr->affordanceCollection();
      else if constexpr (std::tuple_size<AffoCollectionOfAssembler>::value == sizeof...(Affordances))
        return Impl::optionalDecider([&]() { return assemblerPtr->boundToAffordanceCollection(); },
                                     [&]() { return assemblerPtr->affordanceCollection(); }, affordancesArg);
      else {
        if (affordancesArg.has_value())
          return affordancesArg.value();
        else
          DUNE_THROW(Dune::InvalidStateException,
                     "Neither the assembler is bound to an affordance collection nor was a value passed");
      }
    }();

    return op(assemblerPtr, req, affordances, dbcOption,funcIndices);
  }

  // template <typename Assembler>
  // static auto op(Assembler&& as, DBCOption dbcOption) {
  //   auto ex = []() {
  //     DUNE_THROW(Dune::InvalidStateException,
  //                "Assembler has to be bound to a fe requirement and an affordance collection before you can call "
  //                "this method");
  //   };
  //   if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
  //                 traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
  //     if (as->boundToRequirement() and as->boundToAffordanceCollection()) {
  //       return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), dbcOption);
  //     } else {
  //       ex();
  //     }
  //   } else {
  //     if (as->boundToRequirement() and as->boundToAffordanceCollection()) {
  //       return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), dbcOption);
  //     } else {
  //       ex();
  //     }
  //   }
  //   __builtin_unreachable();
  // }

  // template <typename Assembler>
  // static auto op(Assembler&& as) {
  //   auto ex = []() {
  //     DUNE_THROW(Dune::InvalidStateException,
  //                "Assembler has to be bound to a fe requirement to an affordance collection and to an "
  //                "DBCOption before you can call "
  //                "this method");
  //   };
  //   if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
  //                 traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
  //     if (not as->bound())
  //       ex();
  //     return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), as->dBCOption());
  //   } else {
  //     if (not as.bound())
  //       ex();
  //     return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), as->dBCOption());
  //   }
  // }

  // template <typename Assembler, typename... Affordances>
  // static auto op(Assembler&& as, AffordanceCollection<Affordances...> affordances,
  //                DBCOption dbcOption = DBCOption::Full) {
  //   auto ex = []() {
  //     DUNE_THROW(Dune::InvalidStateException,
  //                "Assembler has to be bound to a fe requirement before you can call "
  //                "this method");
  //   };

  // if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
  //               traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
  //   if (not as->boundToRequirement())
  //     ex();
  //   return op(std::forward<Assembler>(as), as->requirement(), affordances, dbcOption);
  // } else {
  //   if (not as.boundToRequirement())
  //     ex();
  //   return op(std::forward<Assembler>(as), as.requirement(), affordances, dbcOption);
  // }
  // }

  // template <typename Assembler>
  // static auto op(Assembler&& as, typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement&
  // req,
  //                DBCOption dbcOption) {
  //   auto ex = []() {
  //     DUNE_THROW(Dune::InvalidStateException,
  //                "Assembler has to be bound to an affordance collection before you can call "
  //                "this method");
  //   };

  // if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
  //               traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
  //   if (not as->boundToAffordanceCollection())
  //     ex();
  //   return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), dbcOption);
  // } else {
  //   if (not as.boundToAffordanceCollection())
  //     ex();
  //   return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), dbcOption);
  // }
  // }
};

} // namespace Ikarus
