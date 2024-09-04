// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <dune/python/pybind11/pybind11.h>
#include <spdlog/fmt/bundled/core.h>
namespace Ikarus::Python {
  void registerSolverInformation(){
      auto solvers = pybind11::module::import("ikarus.solvers");
    namespace py = pybind11;
      py::class_<NonLinearSolverInformation>(solvers, "NonLinearSolverInformation")
        .def(py::init<>())  // Default constructor
        .def_readonly("success", &NonLinearSolverInformation::success)
        .def_readonly("residualNorm", &NonLinearSolverInformation::residualNorm)
        .def_readonly("correctionNorm", &NonLinearSolverInformation::correctionNorm)
        .def_readonly("iterations", &NonLinearSolverInformation::iterations)
        .def("__bool__", [](const NonLinearSolverInformation &self) {
            return static_cast<bool>(self);
        }) .def("__repr__", [](const NonLinearSolverInformation &self) {
            std::string status = self.success ? "Optimization terminated successfully." : "Optimization failed.";
            return fmt::format(
                "{}\n"
                "         Residual Norm: {:.6e}\n"
                "         Correction Norm: {:.6e}\n"
                "         Iterations: {}\n",
                status,
                self.residualNorm,
                self.correctionNorm,
                self.iterations
            );
        });
  }
}