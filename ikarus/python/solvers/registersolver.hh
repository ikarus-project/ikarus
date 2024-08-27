// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/solver/nonlinearsolver/trustregion.hh>


namespace Ikarus::Python {

/**
 * @brief Registers the trust region non-linear solver class in Python.
 *
 * @tparam TR The solver class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class TR, class... options>
void registerTrustRegion(pybind11::handle scope, pybind11::class_<TR, options...> cls) {
  namespace py = pybind11;

  using NonLinearOperator = TR::NonLinearOperator;
  using UpdateFunction    = TR::UpdateFunction;

  cls.def(pybind11::init([](const NonLinearOperator& nonLinearOperator, UpdateFunction updateFunction) {
            return new TR(nonLinearOperator, updateFunction);
          }),
          py::arg("nonLinearOperator"), py::arg("updateFunction") = Ikarus::utils::UpdateDefault{});

  cls.def(
      "setup",
      [](TR& self, py::dict dict) {
    TrustRegionSettings settings;
    settings.populate(dict);
    self.setup(settings);
        // settingsFromDict<int>( "verbosity" );
        // settingsFromDict<double>( "maxtime" );
        // settingsFromDict<int>( "minIter" );
        // settingsFromDict<int>( "maxIter" );
        // settingsFromDict<int>( "debug" );
        // settingsFromDict<double>( "grad_tol" );
        // settingsFromDict<double>( "corr_tol" );
        // settingsFromDict<double>( "rho_prime" );
        // settingsFromDict<bool>( "useRand" );
        // settingsFromDict<double>( "rho_reg" );
        // settingsFromDict<double>( "Delta_bar" );
        // settingsFromDict<double>( "Delta0" );
},
      R"(
        This function sets up the solver with settings provided
        through a Python dictionary.

        Parameters:
            dict (dict): Python dictionary containing trust region settings.

        Notes:
            Supported trust region settings:
            - "verbosity": Verbosity level (int).
            - "maxtime": Maximum allowable time for solving (float).
            - "minIter": Minimum number of outer iterations (int).
            - "maxIter": Maximum number of outer iterations (int).
            - "debug": Debugging flag (int).
            - "grad_tol": Gradient tolerance (float).
            - "corr_tol": Correction tolerance (float).
            - "rho_prime": Rho prime value (float).
            - "useRand": Flag for using random correction predictor (bool).
            - "rho_reg": Regularization value for rho (float).
            - "Delta_bar": Maximum trust region radius (float).
            - "Delta0": Initial trust region radius (float).

        Example:
            ```python
            solver = MySolver()
            solver.setup({'verbosity': 2, 'maxtime': 10.0, 'minIter': 5}));
        )");
}
} // namespace Ikarus::Python