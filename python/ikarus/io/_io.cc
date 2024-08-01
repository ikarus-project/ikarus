// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "../pythonhelpers.hh"

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/io/vtkdatatag.hh>

PYBIND11_MODULE(_io, m) {
  namespace py = pybind11;
  using namespace pybind11::literals;

  using namespace Ikarus::Vtk;
  using namespace Ikarus;
  ENUM_BINDINGS(DataTag);
}