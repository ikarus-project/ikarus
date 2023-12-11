// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <Python.h>
#include <cstddef>
#include <memory>
#include <utility>

#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/basistags.hh>

#include <autodiff/forward/real/real.hpp>

#include <ikarus/utils/basis.hh>

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
