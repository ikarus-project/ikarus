// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/python/pybind11/pybind11.h>

namespace Ikarus::Python {

/**
 * @brief Registers a NonLinearElasticPre class in Python.
 *
 * @tparam NonLinearElasticPre The NonLinearElasticPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class NonLinearElasticPre, class... options>
void registerNonLinearElasticPre(pybind11::handle scope, pybind11::class_<NonLinearElasticPre, options...> cls) {
  using Material = typename NonLinearElasticPre::Material;
  cls.def(pybind11::init([](const Material& mat) { return new NonLinearElasticPre(mat); }));
}

/**
 * @brief Registers a LinearElasticPre class in Python.
 *
 * @tparam LinearElasticPre The LinearElasticPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class LinearElasticPre, class... options>
void registerLinearElasticPre(pybind11::handle scope, pybind11::class_<LinearElasticPre, options...> cls) {
  cls.def(pybind11::init([](double emod, double nu) { return new LinearElasticPre({emod, nu}); }));
}

/**
 * @brief Registers a KirchhoffLoveShellPre class in Python.
 *
 * @tparam KirchhoffLoveShellPre The KirchhoffLoveShellPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class KirchhoffLoveShellPre, class... options>
void registerKirchhoffLoveShellPre(pybind11::handle scope, pybind11::class_<KirchhoffLoveShellPre, options...> cls) {
  cls.def(pybind11::init(
      [](const double& E, const double& nu, const double& h) { return new KirchhoffLoveShellPre({E, nu}, h); }));
}

/**
 * @brief Registers an EnhancedAssumedStrainsPre class in Python.
 *
 * @tparam EASPre The EnhancedAssumedStrainsPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class EASPre, class... options>
void registerEnhancedAssumedStrainsPre(pybind11::handle scope, pybind11::class_<EASPre, options...> cls) {
  cls.def(pybind11::init([](int numberOfParameter) { return new EASPre(numberOfParameter); }));
}

/**
 * @brief Registers a NeumannBoundaryLoadPre class in Python.
 *
 * @tparam NeumannBoundaryLoadPre The NeumannBoundaryLoadPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
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

/**
 * @brief Registers a VolumeLoadPre class in Python.
 *
 * @tparam VolumeLoadPre The VolumeLoadPre class.
 * @tparam options Additional options for the pybind11 class.
 *
 * @param scope Python handle to the module or class scope.
 * @param cls The pybind11 class to register.
 */
template <class VolumeLoadPre, class... options>
void registerVolumeLoadPre(pybind11::handle scope, pybind11::class_<VolumeLoadPre, options...> cls) {
  using LoadFunction = std::function<Eigen::Vector<double, VolumeLoadPre::worldDim>(
      Dune::FieldVector<double, VolumeLoadPre::worldDim>, const double&)>;
  cls.def(pybind11::init([](LoadFunction volumeLoad) { return new VolumeLoadPre(volumeLoad); }));
}

} // namespace Ikarus::Python
