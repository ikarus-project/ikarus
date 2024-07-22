// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <dune/common/classname.hh>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/nonlinopfactory.hh>
#include <ikarus/utils/integersequence.hh>
namespace Ikarus::Python {

template <class NLO, class... options, typename T,size_t sizeOfIndices>
void registerNonLinearOperator(pybind11::handle scope, pybind11::class_<NLO, options...> cls, std::array<T,sizeOfIndices> baseIndices) ;

namespace Impl {
  template <typename NLO, int... i>
  struct SubOperatorHelper
  {
    using SubOperator = decltype(std::declval<NLO>().template subOperator<i...>());
  };

  template<typename T,T... I>
  std::string indexSeqToString(std::integer_sequence<T,I...>){
    return (std::to_string(I) + ...+"");
  }

   template<typename NLO>
   auto indexSeqOfFunctionIndices = Dune::makeArrayFromIndexSequence(std::make_index_sequence<std::tuple_size_v<typename NLO::FunctionTuple>>{});



  template <typename NLO,typename NLOBASE=NLO, typename T,size_t sizeOfIndices,size_t sizeOfSubIndices>
  auto maybeRegisterNonlinearOperator(pybind11::handle scope, std::array<T,sizeOfIndices> baseIndices, std::array<T,sizeOfSubIndices> subIndices) {
    std::cout<<"====================="<<std::endl;
    std::cout<<"maybeRegisterNonlinearOperator"<<std::endl;
    std::cout<<"NLO: "<<Dune::className<NLO>()<<std::endl;
    auto includes = Dune::Python::IncludeFiles{"ikarus/utils/nonlinearoperator.hh"};

    /* The names of the nonlinear operators are generated by concatenating the names of the subfunctions and the indices of the subfunctions.
    This is needed to not artificially insert the same nonlinear operator multiple times with different names (inefficient)
    and even worse not try to insert different nonlinear operators with the same name (Python registry import error).
    We do basically here
    NonlinearOperator_012 -> NonlinearOperator_012
NonlinearOperator_012_0 -> NonlinearOperator_0
NonlinearOperator_012_01 -> NonlinearOperator_01
NonlinearOperator_012_01_0 -> NonlinearOperator_0
NonlinearOperator_012_01_1 -> NonlinearOperator_1
NonlinearOperator_012_1 -> NonlinearOperator_1
NonlinearOperator_012_12 -> NonlinearOperator_12
NonlinearOperator_012_12_0  -> NonlinearOperator_1
NonlinearOperator_012_12_1 -> NonlinearOperator_2
NonlinearOperator_012_2  -> NonlinearOperator_2
        */

     std::string name = "NonLinearOperator_";
     std::cout<<"baseIndices"<<std::endl;
     for(auto i: baseIndices)
     std::cout<<" "<<i<<std::endl;

             std::cout<<"subIndices"<<std::endl;
     for(auto i: subIndices)
     std::cout<<" "<<i<<std::endl;
     std::array<T,std::tuple_size_v<typename NLO::FunctionTuple>> subFunctions;
    for(auto i=0;i<subIndices.size();i++)
        {
            subFunctions[i]=baseIndices[subIndices[i]];
            name+= std::to_string(subFunctions[i]);
            }
    auto [clsNonLinearOperator, notRegisteredNonLinOp] =
        Dune::Python::insertClass<NLO>(scope, name, pybind11::dynamic_attr(),
                                       Dune::Python::GenerateTypeName(Dune::className<NLO>()), includes);

    if (notRegisteredNonLinOp)
      {
        std::cout<<"register "+name<<std::endl;
        registerNonLinearOperator(scope, clsNonLinearOperator,subFunctions);
        }else
         std::cout<<"Already registered NonLinearOperator "+name<<std::endl;
    std::cout<<"====================="<<std::endl;
    return clsNonLinearOperator;
  }

    template <typename NLO,typename NLOBASE=NLO, typename T,size_t sizeOfIndices>
  auto maybeRegisterNonlinearOperator(pybind11::handle scope, std::array<T,sizeOfIndices> baseIndices) {
    return maybeRegisterNonlinearOperator<NLO,NLOBASE>(scope,baseIndices,baseIndices);
  }

  template <int... i, typename NLO, class... options, typename T,size_t sizeOfIndices>
  auto maybeRegisterNonlinearSubOperator(pybind11::handle scope, pybind11::class_<NLO, options...> cls, std::array<T,sizeOfIndices> baseIndices) {
    std::cout<<"====================="<<std::endl;
    std::cout<<"maybeRegisterNonlinearSubOperator"<<std::endl;
    std::cout<<"NLO: "<<Dune::className<NLO>()<<std::endl;
    std::cout<<"NLOSUB: "<<Dune::className<typename SubOperatorHelper<NLO, i...>::SubOperator>()<<std::endl;
    std::cout<<"i: "<<(std::to_string(i) + ...)<<std::endl;
    std::array subIndices = Dune::makeArrayFromIndexSequence<T,i...>();
    auto clsSub = maybeRegisterNonlinearOperator<typename SubOperatorHelper<NLO, i...>::SubOperator,NLO>(
        scope,baseIndices,subIndices);
    cls.def(("__subOperator" + (std::to_string(i) + ...)).c_str(),
            [](NLO& self) { return self.template subOperator<i...>(); });
    std::cout<<"====================="<<std::endl;

    return clsSub;
  }
} // namespace Impl

/**
 * @brief Registers a nonlinear operator class in Python.
 *
 * @tparam NLO The nonlinear operator class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class NLO, class... options, typename T,size_t sizeOfIndices>
void registerNonLinearOperator(pybind11::handle scope, pybind11::class_<NLO, options...> cls, std::array<T,sizeOfIndices> baseIndices) {
  namespace py = pybind11;

  using ParameterValues = NLO::ParameterValues;
  using FunctionTuple   = NLO::FunctionTuple;

  cls.def(pybind11::init([](const FunctionTuple& fs, const ParameterValues& parameters) {
    auto lambdaParameter = [](auto&&... v) -> decltype(auto) { return parameter(v...); };
    auto lambdaFunc      = [](auto&&... v) -> decltype(auto) { return functions(v...); };
    return new NLO(std::apply(lambdaFunc, fs), std::apply(lambdaParameter, parameters));
  }));

  cls.def("updateAll", &NLO::updateAll, "Update all functions values");

  constexpr int numberOfFunctions  = std::tuple_size_v<FunctionTuple>;
  constexpr int numberOfParameters = std::tuple_size_v<ParameterValues>;

  cls.def_property_readonly_static("numberOfFunctions", [](py::object) { return numberOfFunctions; });

  cls.def_property_readonly_static("numberOfParameters", [](py::object) { return numberOfParameters; });
  cls.def(
      "update",
      [numberOfFunctions](NLO& self, size_t index) {
        if (index >= numberOfFunctions)
          throw py::index_error();
        Dune::Hybrid::switchCases(Dune::Hybrid::integralRange(Dune::index_constant<numberOfFunctions>()), index,
                                  [&](auto i) { self.template update<i>(); });
      },
      py::arg("index"), "Update the function value of the given index");

  if constexpr (numberOfFunctions > 0) {
    cls.def("value", &NLO::value, py::return_value_policy::reference_internal);
    cls.def_property_readonly_static(
        "valueCppTypeName", [](py::object) { return Dune::className<typename NLO::template FunctionReturnType<0>>(); });
  }


  if constexpr (numberOfFunctions > 1) {
    cls.def("derivative", &NLO::derivative, py::return_value_policy::reference_internal);
    cls.def_property_readonly_static("derivativeCppTypeName", [](py::object) {
      return Dune::className<typename NLO::template FunctionReturnType<1>>();
    });
    Impl::maybeRegisterNonlinearSubOperator<0>(scope, cls,baseIndices);
    Impl::maybeRegisterNonlinearSubOperator<1>(scope, cls,baseIndices);
  }

  if constexpr (numberOfFunctions > 2) {
    cls.def("secondDerivative", &NLO::secondDerivative, py::return_value_policy::reference_internal);
    cls.def_property_readonly_static("secondDerivativeCppTypeName", [](py::object) {
      return Dune::className<typename NLO::template FunctionReturnType<2>>();
    });
    Impl::maybeRegisterNonlinearSubOperator<0, 1>(scope, cls,baseIndices);
    Impl::maybeRegisterNonlinearSubOperator<1, 2>(scope, cls,baseIndices);
    Impl::maybeRegisterNonlinearSubOperator<2>(scope, cls,baseIndices);
  }

  if constexpr (numberOfFunctions > 0) {
    cls.def("lastParameter", &NLO::lastParameter, py::return_value_policy::reference_internal);
    cls.def_property_readonly_static("lastParameterCppTypeName", [](py::object) {
      return Dune::className<typename NLO::template ParameterValue<numberOfParameters - 1>>();
    });
    cls.def("firstParameter", &NLO::firstParameter, py::return_value_policy::reference_internal);
    cls.def_property_readonly_static("firstParameterCppTypeName", [](py::object) {
      return Dune::className<typename NLO::template ParameterValue<0>>();
    });
  }
  if constexpr (numberOfFunctions > 1) {
    cls.def("secondParameter", &NLO::secondParameter, py::return_value_policy::reference_internal);
    cls.def_property_readonly_static("secondParameterCppTypeName", [](py::object) {
      return Dune::className<typename NLO::template ParameterValue<1>>();
    });
  }
}

template <class NLO, class... options>
void registerNonLinearOperator(pybind11::handle scope, pybind11::class_<NLO, options...> cls) {
  return registerNonLinearOperator(scope, cls, Ikarus::Python::Impl::indexSeqOfFunctionIndices<NLO>);
}

template <typename A>
struct NonLinearOperatorFactoryWrapper
{
  using Assembler = A;
  Assembler as;
};

/**
 * @brief Registers the nonlinear operator factory in Python.
 *
 * @tparam NLOFW The factory wrapper struct from above.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class NLOFW, class... options>
void registerNonLinearOperatorFactory(pybind11::handle scope, pybind11::class_<NLOFW, options...> cls) {
  namespace py = pybind11;

  using Assembler            = NLOFW::Assembler;
  using AffordanceCollection = AffordanceCollection<ScalarAffordance, VectorAffordance, MatrixAffordance>;

  cls.def(pybind11::init([](const Assembler& as) { return new NLOFW{as}; }));

  cls.def(pybind11::init([](const Assembler& as) { return new NLOFW{.as = as}; }));



  auto registerNonLinearOperatorL = [&scope]<size_t... i>(std::index_sequence <i ...> iSeq) {

 Impl::maybeRegisterNonlinearOperator<decltype(NonLinearOperatorFactory::template op<i...>(std::declval<Assembler&>()))>(scope,Dune::makeArrayFromIndexSequence(iSeq));
  };

  registerNonLinearOperatorL(std::index_sequence<0,1,2>());



  // cls.def(
  //     "op",
  //     [](NLOFW& self,
  //        std::optional<
  //            std::reference_wrapper<typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement>>
  //            req,
  //        std::optional<AffordanceCollection> affordances,
  //        std::optional<DBCOption> enforce) { return NonLinearOperatorFactory::op(self.as, req, affordances, enforce); },
  //     "Create a nonlinear operator from the given assembler and the optional fe requirements, affordances and "
  //     "enforcing dbc option.",
  //     py::arg("requirement") = py::none(), py::arg("affordances") = py::none(),
  //     py::arg("enforcingDBCOption") = py::none());



}

} // namespace Ikarus::Python