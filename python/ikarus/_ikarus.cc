// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>


#include <ikarus/python/ikarus/enums.hh>
#include <ikarus/affordances.hh>

#include <dune/python/pybind11/pybind11.h>
int add(int i, int j) {
  return i + j;
}
PYBIND11_MODULE( _ikarus, m )
{
  // enumeration types from dune-grid
  m.def("add", &add, "A function which adds two numbers",
        pybind11::arg("i"), pybind11::arg("j"));

  pybind11::enum_< Ikarus::Bla > scalarAffordances( m, "Bla" );
  scalarAffordances.value( "foo", Ikarus::Bla::foo );
  scalarAffordances.value( "bar", Ikarus::Bla::bar );
}
