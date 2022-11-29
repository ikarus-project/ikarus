/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#pragma once
#include <Python.h>

#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

#include <autodiff/forward/real/real.hpp>
namespace Ikarus {
  template <typename... Args>
  auto makeSharedBasis(Args&&... args) {
    using namespace Dune::Functions::BasisFactory;
    using DuneBasis = decltype(makeBasis(std::forward<Args>(args)...));
    return std::make_shared<DuneBasis>(makeBasis(std::forward<Args>(args)...));
  }

  template <typename... Args>
  auto makeConstSharedBasis(Args&&... args) {
    using namespace Dune::Functions::BasisFactory;
    using DuneBasis = decltype(makeBasis(std::forward<Args>(args)...));
    return std::make_shared<const DuneBasis>(makeBasis(std::forward<Args>(args)...));
  }
}  // namespace Ikarus

namespace Python {
  // *****************************************************************************
  // specializations of Conversion that do the PyObject* <-> C++-type conversion
  // *****************************************************************************

  // conversion of autodiff::Real
  template <std::size_t order, class T>
  struct Conversion<autodiff::Real<order, T>> {
    enum { useDefaultConstructorConversion = true };
    static void toC(PyObject* list, autodiff::Real<order, T>& v) {
      auto rlist = Reference(Imp::inc(list));
      // Reference is needed to enable the ".get()" function and "Imp::inc" is
      // needed since Reference owns the PyObject pointer and decrements it at the end of the scope
      // Imp::inc artificially increases the reference counter by one.
      // When we return from this function, these two cancel out and the PyObject* is as before

      auto wF = Callable(rlist.get("__getitem__"));
      for (std::size_t i = 0; i < order + 1; ++i)
        v[i] = PyFloat_AsDouble(wF(i));
    }

    static PyObject* toPy(const autodiff::Real<order, T>& v) {
      auto pyMain           = Python::main();
      Python::Module module = pyMain.import("autodiff");

      auto real1stClass = module.get("real1st");
      auto real1st      = Callable(Imp::inc(real1stClass))();
      auto wF           = Callable(Imp::inc(real1st).get("__setitem__"));
      for (std::size_t i = 0; i < order + 1; ++i)
        wF(i, v[i]);  // real1st.__setitem__(i,v[i]);

      return real1st;
    }
  };
}  // namespace Python