// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file linearelastic.hh
 * \brief Python bindings for the linear elastic element
 */

#pragma once

#include "registerelement.hh"

#include <dune/fufem/boundarypatch.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/utils/basis.hh>

namespace Ikarus::Python {

/**
 * \brief Register Python bindings for a Linear Elastic class.
 *
 * This function registers Python bindings for a LinearElastic class, allowing it to be used in Python scripts.
 * The registered class will have several initializers with different sets of parameters.
 *
 * \tparam LinearElastic The LinearElastic class to be registered.
 * \tparam options Variadic template parameters for additional options when defining the Python class.
 *
 * \param scope A Pybind11 handle representing the Python scope where the class should be registered.
 * \param cls The Pybind11 class template to be used for registering the LinearElastic class.
 *
 * \ingroup pythonbindings
 */
template <class LinearElastic, class... options>
void registerLinearElastic(pybind11::handle scope, pybind11::class_<LinearElastic, options...> cls) {
  registerElement(scope, cls);
}

} // namespace Ikarus::Python
