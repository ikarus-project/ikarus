// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/eigen.h>

#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>


namespace Ikarus::Python {

/**
 * @brief Registers the Newton Raphson non-linear solver class in Python.
 *
 * @tparam NR The solver class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class NR, class... options>
void registerNewtonRaphson(pybind11::handle scope, pybind11::class_<NR, options...> cls) {
  namespace py = pybind11;

  using NonLinearOperator = NR::NonLinearOperator;
  using UpdateFunction    = NR::UpdateFunction;
  using LinearSolver = NR::LinearSolver;

  cls.def(pybind11::init([](const NonLinearOperator& nonLinearOperator,Ikarus::SolverTypeTag solverTag, UpdateFunction updateFunction) {
            return new NR(nonLinearOperator, Ikarus::LinearSolver(solverTag), std::move(updateFunction));
          }),
          py::arg("nonLinearOperator"), py::arg("linearSolverTag"), py::arg("updateFunction") = UpdateFunction(Ikarus::utils::UpdateDefault{}));

  cls.def(
      "setup",
      [](NR& self, py::dict dict) {
    NewtonRaphsonSettings settings;
    settings.populate(dict);
    self.setup(settings);
},
      R"(
        This function sets up the solver with settings provided
        through a Python dictionary.

        Parameters:
            dict (dict): Python dictionary containing Newton Raphson settings.

        Notes:
            Supported Newton Raphson settings:
            - "maxIter": Maximum number of outer iterations (int).
            - "grad_tol": Gradient tolerance (float).

        Example:
            ```python
            solver = NewtonRaphson()
            solver.setup({'maxIter': 100, 'grad_tol': 1e-8}));
        )");

          cls.def(
      "solve",
      [](NR& self, const typename NR::CorrectionType& dx) {
          self.solve(dx);
      }, py::arg("predictor") );

  cls.def(
      "solve",
      [](NR& self) {
          self.solve();
      });

      cls.def("nonLinearOperator", &NR::nonLinearOperator, py::return_value_policy::reference_internal);
}


} // namespace Ikarus::Python