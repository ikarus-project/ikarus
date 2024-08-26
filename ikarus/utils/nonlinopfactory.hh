// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinopfactory.hh
 * \brief Contains the generic NonLinearOperatorFactory class.
 */

#pragma once

#include <utility>

#include <dune/common/indices.hh>
#include <dune/common/integersequence.hh>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/traits.hh>

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
  struct DummyEmpty
  {
  };
  template <int index, typename Assembler, typename... Affordances>
  requires(index < 3 and index >= 0)
  static auto function(std::shared_ptr<Assembler> assemblerPtr, AffordanceCollection<Affordances...> affordancesInput,
                       DBCOption dbcOption) {
    using FERequirement = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;
    struct DummyEmpty
    {
    };

    auto affordances = [&]() {
      if constexpr (sizeof...(Affordances) == 0)
        return assemblerPtr->affordanceCollection();
      else
        return affordancesInput;
    }();

    // Since it is not possible to have [[no_unique_address]] with a lambda, we have to use a dummy type, to remove the
    // lambda overhead capturing dbcOption for the scalar function
    struct DummyLambda
    {
      decltype(auto) operator()(typename FERequirement::SolutionVectorType& globalSol,
                                typename FERequirement::ParameterType& parameter) const {
        FERequirement req;
        req.insertGlobalSolution(globalSol).insertParameter(parameter);

        if constexpr (index == 0)
          return assembler->scalar(req, affordancesArg.scalarAffordance());
        else if constexpr (index == 1)
          return assembler->vector(req, affordancesArg.vectorAffordance(), dbcOption);
        else if constexpr (index == 2)
          return assembler->matrix(req, affordancesArg.matrixAffordance(), dbcOption);
      }

      std::shared_ptr<Assembler> assembler;
      std::remove_reference_t<decltype(affordances)> affordancesArg;
      [[no_unique_address]] std::conditional_t<index == 0, DummyEmpty, DBCOption> dbcOption;
    };
    DummyLambda result;
    result.assembler      = assemblerPtr;
    result.affordancesArg = affordances;
    if constexpr (index != 0)
      result.dbcOption = dbcOption;

    return result;
  }

  template <size_t... funcs, typename Assembler, typename... Affordances>
  static auto op_impl2(Assembler&& as,
                       typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement& req,
                       AffordanceCollection<Affordances...> affordances, DBCOption dbcOption,
                       std::index_sequence<funcs...> funcIndices = {}) {}

  template <typename Assembler, typename... Args>
  static auto op(Assembler&& as, Args&&... args) {
    // Define the default values
    using ReqType = std::optional<
        std::reference_wrapper<typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement>>;
    ReqType reqArg                                       = std::nullopt;
    std::optional<AffordanceCollection<>> affordancesArg = std::nullopt;
    std::optional<DBCOption> dbCOptionArg                = std::nullopt;
    std::index_sequence<> funcIndicesDefault             = {};

    auto aCCorrectIndex          = Dune::index_constant<1>{};
    auto funcIndicesCorrectIndex = Dune::index_constant<3>{};
    auto argumentDefaultTuple    = std::make_tuple(reqArg, affordancesArg, dbCOptionArg, funcIndicesDefault);
    auto argumentTupleRaw        = std::make_tuple(std::forward<Args>(args)...);

    constexpr auto AIndex = Ikarus::utils::find_if(argumentTupleRaw, [](const auto& arg) {
      return Ikarus::IsAffordanceCollection_v<std::remove_cvref_t<decltype(arg)>>;
    });

    constexpr auto FIIndex   = Ikarus::utils::find_if(argumentTupleRaw, [](const auto& arg) {
      return Ikarus::traits::IsIntegerSequence_v<std::remove_cvref_t<decltype(arg)>>;
    });
    auto argumentTupleLambda = [&](auto tuple, auto tupleRaw, auto correctIndex, auto index,
                                   auto wrapInOptional) constexpr -> decltype(auto) {
      if constexpr (index < sizeof...(Args)) {
        using NewType = typename std::tuple_element<index, decltype(tupleRaw)>::type;
        auto result =
            Ikarus::traits::ReplaceTypeAtPos_t<decltype(tuple), correctIndex,
                                               std::conditional_t<wrapInOptional, std::optional<NewType>, NewType>>{};
        Dune::Hybrid::forEach(std::make_index_sequence<std::tuple_size_v<decltype(result)>>{}, [&](auto i) {
          if constexpr (i != correctIndex)
            std::get<i>(result) = std::get<i>(tuple);
          else
            std::get<correctIndex>(result) = std::get<index>(tupleRaw);
        });

        return result;
      } else
        return tuple;
    };
    auto argumentTupleT = argumentTupleLambda(argumentDefaultTuple, argumentTupleRaw, aCCorrectIndex,
                                              Dune::index_constant<AIndex>{}, std::bool_constant<true>{});
    // static_assert(std::tuple_size_v<decltype()> ==4);
    auto argumentTuple = argumentTupleLambda(argumentTupleT, argumentTupleRaw, funcIndicesCorrectIndex,
                                             Dune::index_constant<FIIndex>{}, std::bool_constant<false>{});

    auto range = std::make_index_sequence<sizeof...(Args)>{};
    Dune::Hybrid::forEach(range, [&](auto i) {
      auto arg  = std::get<i>(argumentTupleRaw);
      using Arg = std::remove_cvref_t<decltype(arg)>;
      std::cout << "Arg:" << Dune::className<Arg>() << std::endl;
      if constexpr (std::is_same_v<Arg,
                                   typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement&>) {
        std::get<0>(argumentTuple) = std::get<i>(argumentTupleRaw);
      } else if constexpr (std::is_same_v<Arg, DBCOption>) {
        std::get<2>(argumentTuple) = std::get<i>(argumentTupleRaw);
      }
    });

    return std::apply(
        [&]<typename... Args2>(Args2&&... args2) {
          return op_impl(std::forward<Assembler>(as), std::forward<Args2>(args2)...);
        },
        argumentTuple);
  }

private:
  template <size_t... funcs, typename Assembler, typename... Affordances>
  static auto op_impl(
      Assembler&& as,
      std::optional<
          std::reference_wrapper<typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement>>
          reqArg                                                         = std::nullopt,
      std::optional<AffordanceCollection<Affordances...>> affordancesArg = std::nullopt,
      std::optional<DBCOption> dbCOptionArg = std::nullopt, std::index_sequence<funcs...> funcIndices = {}) {
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

    constexpr int funcIndexSize = funcIndices.size();
    static_assert(Dune::equal(Dune::sorted(funcIndices), funcIndices),
                  "The function indices you request have to be sorted");
    static_assert(funcIndexSize < 4, "The number of function indices you request have to be less than 4");
    static_assert(Dune::filter([](auto i) { return i < 3; }, funcIndices).size() == funcIndexSize,
                  "The function indices you request have to be less than 3");
    using namespace Dune::Indices;

    constexpr bool provideScalar =
        (affordances.hasScalarAffordance and (funcIndexSize == 0 or Dune::contains(funcIndices, _0)));
    constexpr bool provideVector =
        (affordances.hasVectorAffordance and (funcIndexSize == 0 or Dune::contains(funcIndices, _1)));
    constexpr bool provideMatrix =
        (affordances.hasMatrixAffordance and (funcIndexSize == 0 or Dune::contains(funcIndices, _2)));

    constexpr std::array<bool, 3> provide{provideScalar, provideVector, provideMatrix};
    auto funcs2 =
        Dune::filter([&](auto i) constexpr { return std::bool_constant<provide[i]>{}; }, std::make_index_sequence<3>{});
    static_assert(funcs2.size() == provideScalar + provideVector + provideMatrix);

    assert(req.populated() && " Before you calls this method you have to pass populated fe requirements");
    auto createNonLinearOp = [&]<size_t... funcs3>(std::index_sequence<funcs3...>) {
      return NonLinearOperator(functions(function<funcs3>(assemblerPtr, affordances, dbcOption)...),
                               parameter(req.globalSolution(), req.parameter()));
    };
    return createNonLinearOp(funcs2);
  }
};

} // namespace Ikarus
