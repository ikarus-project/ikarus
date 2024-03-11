// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file traction.hh
 * \brief Python bindings for the traction pre
 */

#pragma once

#include <dune/python/common/typeregistry.hh>
#include <dune/python/functions/globalbasis.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/mechanics/loads/traction.hh>

namespace Ikarus::Python {

/**
 * \brief Register Python bindings for the NeumannBoundaryLoadPre class.
 *
 * This function registers Python bindings for a NeumannBoundaryLoadPre class, allowing it to be used in Python scripts.
 *
 * \tparam NeumannBoundaryLoadPre The NeumannBoundaryLoadPre class to be registered.
 * \tparam options Variadic template parameters for additional options when defining the Python class.
 *
 * \param scope A Pybind11 handle representing the Python scope where the class should be registered.
 * \param cls The Pybind11 class template to be used for registering the KirchhoffLoveShell class.
 *
 * \ingroup pythonbindings
 */
template <class NeumannBoundaryLoadPre, class... options>
void registerNeumannBoundaryLoadPre(pybind11::handle scope, pybind11::class_<NeumannBoundaryLoadPre, options...> cls) {
  using BoundaryPatchType = typename NeumannBoundaryLoadPre::BoundaryPatchType;
  using GridView          = typename NeumannBoundaryLoadPre::GridView;

  using LoadFunction = std::function<Eigen::Vector<double, NeumannBoundaryLoadPre::worldDim>(
      Dune::FieldVector<double, NeumannBoundaryLoadPre::worldDim>, const double&)>;
  cls.def(pybind11::init([](const BoundaryPatchType& patch, LoadFunction volumeLoad) {
            return new NeumannBoundaryLoadPre(&patch, volumeLoad);
          }),
          pybind11::keep_alive<1, 2>());
}

} // namespace Ikarus::Python
