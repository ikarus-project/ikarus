
// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file pathfollowingfunctionspythonwrapper.hh
 * \brief Defines structures and methods related to subsidiary functions for control routines.
 *
 */

 #include <ikarus/controlroutines/pathfollowingfunctions.hh>
 #include <ikarus/utils/nonlinearoperator.hh>
 #include <ikarus/controlroutines/pathfollowingfunctionsinterface.hh>

 namespace Ikarus::Python{

template<NLO>
struct PySubsidaryFunction :  SubsidaryFunction<NLO> {
   using SubsidaryFunction<NLO>::SubsidaryFunction<NLO>;

    void operator()(SubsidiaryArgs& args) {
    PYBIND11_OVERRIDE_PURE_NAME(void, SubsidaryFunction<NLO>, "__call__",operator(), args);
         }
    void initialPrediction(NLO& nonLinearOperator, SubsidiaryArgs& args) {
         PYBIND11_OVERRIDE_PURE(void, SubsidaryFunction<NLO>,initialPrediction, nonLinearOperator,args);
    }
    void intermediatePrediction(NLO& nonLinearOperator, SubsidiaryArgs& args) {
                 PYBIND11_OVERRIDE_PURE(void, SubsidaryFunction<NLO>,intermediatePrediction, nonLinearOperator,args);
         }
    std::string name() const {
                        PYBIND11_OVERRIDE_PURE(std::string, SubsidaryFunction<NLO>,name);
        }
};

template <class NLO>
void registerSubsidaryFunction(pybind11::handle scope) {
    using SF = SubsidaryFunction<NLO>;
    auto includes              = Dune::Python::IncludeFiles{"ikarus/controlroutines/pathfollowingfunctionsinterface.hh"};
    auto [lv, isNotRegistered] = Dune::Python::insertClass<SF,PySubsidaryFunction<NLO>>(
        scope, "SubsidaryFunction", Dune::Python::GenerateTypeName(Dune::className<SF>()), includes);
    if (isNotRegistered) {
        lv.def(pybind11::init<SF>())
        .def("__call__", &SF::operator())
        .def("initialPrediction", &SF::initialPrediction)
        .def("intermediatePrediction", &SF::intermediatePrediction)
        .def("name", &SF::name);
    }
}

} // namespace Python
