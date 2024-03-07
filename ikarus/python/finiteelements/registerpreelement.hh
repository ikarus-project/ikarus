// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/python/pybind11/pybind11.h>

namespace Ikarus::Python {

template <class NonLinearElasticPre, class... options>
void registerNonLinearElasticPre(pybind11::handle scope, pybind11::class_<NonLinearElasticPre, options...> cls) {
  using Material = NonLinearElasticPre::Material;
  cls.def(pybind11::init([](const Material& mat) { return new NonLinearElasticPre(mat); }));
}

template <class LinearElasticPre, class... options>
void registerLinearElasticPre(pybind11::handle scope, pybind11::class_<LinearElasticPre, options...> cls) {
  cls.def(pybind11::init([](double emod, double nu) { return new LinearElasticPre({emod, nu}); }));
}

template <class KirchhoffLoveShellPre, class... options>
void registerKirchhoffLoveShellPre(pybind11::handle scope, pybind11::class_<KirchhoffLoveShellPre, options...> cls) {
  cls.def(pybind11::init(
      [](const double& E, const double& nu, const double& h) { return new KirchhoffLoveShellPre({E, nu}, h); }));
}

template <class EASPre, class... options>
void registerEnhancedAssumedStrainsPre(pybind11::handle scope, pybind11::class_<EASPre, options...> cls) {
  cls.def(pybind11::init([](int numberOfParameter) { return new EASPre(numberOfParameter); }));
}

} // namespace Ikarus::Python
