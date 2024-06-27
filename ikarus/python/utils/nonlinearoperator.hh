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
namespace Ikarus::Python {

template <class NLO, class... options>
void registerNonLinearOperator(pybind11::handle scope, pybind11::class_<NLO, options...> cls) ;

namespace Impl
{
  template<typename NLO,int... i>
  struct SubOperatorHelper
  {
    using SubOperator = decltype(std::declval<NLO>().template subOperator<i...>());
  };

  template<typename NLO, class... Args >
  auto maybeRegisterNonlinearOperator(pybind11::handle scope, std::string name,Args... args )
  {

      auto includes = Dune::Python::IncludeFiles{"ikarus/utils/nonlinearoperator.hh"};
    auto [clsNonLinearOperator, notRegisteredNonLinOp] = Dune::Python::insertClass<NLO>(
    scope,name , pybind11::dynamic_attr(),args..., Dune::Python::GenerateTypeName(Dune::className<NLO>()), includes);

    if (notRegisteredNonLinOp)
      registerNonLinearOperator(scope,clsNonLinearOperator);
      return clsNonLinearOperator;
  }

  template<int... i,typename NLO, class... options>
  auto maybeRegisterNonlinearSubOperator(pybind11::handle scope, pybind11::class_<NLO, options...> cls)
  {
  auto clsSub= maybeRegisterNonlinearOperator<typename SubOperatorHelper<NLO,i...>::SubOperator>(cls,("NonlinearOperator"+(std::to_string(i)+...)));
  cls.def(("__subOperator"+(std::to_string(i)+...)).c_str(), []( NLO& self) {return self.template subOperator<i...>();});
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
template <class NLO, class... options>
void registerNonLinearOperator(pybind11::handle scope, pybind11::class_<NLO, options...> cls) {
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

  cls.def_property_readonly_static("numberOfFunctions", [](py::object) {
     return numberOfFunctions;
   });

     cls.def_property_readonly_static("numberOfParameters", [](py::object) {
     return numberOfParameters;
   });
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
     Impl::maybeRegisterNonlinearSubOperator<0>(scope,cls);
  }
  if constexpr (numberOfFunctions > 1) {
    cls.def("derivative", &NLO::derivative, py::return_value_policy::reference_internal);
    cls.def_property_readonly_static("derivativeCppTypeName", [](py::object) {
      return Dune::className<typename NLO::template FunctionReturnType<1>>();
    });
    Impl::maybeRegisterNonlinearSubOperator<1>(scope,cls);
  }

  if constexpr (numberOfFunctions > 2) {
    cls.def("secondDerivative", &NLO::secondDerivative, py::return_value_policy::reference_internal);
    cls.def_property_readonly_static("secondDerivativeCppTypeName", [](py::object) {
      return Dune::className<typename NLO::template FunctionReturnType<2>>();
    });
      Impl::maybeRegisterNonlinearSubOperator<0,1>(scope,cls);
      Impl::maybeRegisterNonlinearSubOperator<1,2>(scope,cls);
      Impl::maybeRegisterNonlinearSubOperator<2>(scope,cls);
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

  using Assembler       = NLOFW::Assembler;
  using AffordanceCollection = AffordanceCollection<ScalarAffordance,VectorAffordance,MatrixAffordance>;

  cls.def(pybind11::init([](const Assembler& as) {
    return new NLOFW{as};
  }));

  cls.def(pybind11::init(
      [](const Assembler& as) { return new NLOFW{.as = as}; }));

  using NonLinearOperator =
      decltype(NonLinearOperatorFactory::op(std::declval<Assembler&>()));

  Impl::maybeRegisterNonlinearOperator<NonLinearOperator>(cls,"NonLinearOperator");

  cls.def(
      "op",
      [](NLOFW& self,
         std::optional<
             std::reference_wrapper<typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement>>
             req,
         std::optional<AffordanceCollection> affordances, std::optional<DBCOption> enforce) {
          return NonLinearOperatorFactory::op(self.as,req,affordances,enforce);
        }
      ,
      "Create a nonlinear operator from the given assembler and the optional fe requirements, affordances and "
      "enforcing dbc option.",
      py::arg("requirement") = py::none(), py::arg("affordances") = py::none(),
      py::arg("enforcingDBCOption") = py::none());
}

} // namespace Ikarus::Python