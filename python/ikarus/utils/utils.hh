// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <ikarus/utils/tensorutils.hh>

void addBindingsToUtils() {
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;

  auto utils = pybind11::module::import("ikarus.utils");

  /**
   * \brief Converts a square 1x1, 2x2 or 3x3 matrix to a Voigt notation vector.
   *  \ingroup tensor
   * \param E Input matrix of size (size x size).
   * \param isStrain Flag indicating whether the conversion is for strain (true) or not (false) (default is true)..
   * \return  Vector with components in Voigt notation vector.
   *
   * This function converts a square matrix to a Voigt notation vector, which contains the unique components of
   * the input matrix.
   *
   * The optional isStrain parameter allows the user to specify whether the conversion is intended for strain
   * calculations. If isStrain is true, the off-diagonal components are multiplied by 2, providing the correct Voigt
   * notation for symmetric strain tensors.
   */
  utils.def(
      "toVoigt",
      [](Eigen::MatrixXd mat, bool isStrain = true) {
        auto callToVoigt = []<typename T>(const T& m, bool isStrain) {
          auto result = toVoigt(m, isStrain);
          return Eigen::MatrixXd(result);
        };

        if (mat.rows() == 1 && mat.cols() == 1) {
          Eigen::Matrix<double, 1, 1> fixedMat = mat;
          return callToVoigt(fixedMat, isStrain);
        } else if (mat.rows() == 2 && mat.cols() == 2) {
          Eigen::Matrix<double, 2, 2> fixedMat = mat;
          return callToVoigt(fixedMat, isStrain);

        } else if (mat.rows() == 3 && mat.cols() == 3) {
          Eigen::Matrix<double, 3, 3> fixedMat = mat;
          return callToVoigt(fixedMat, isStrain);
        } else {
          DUNE_THROW(Dune::IOError, "toVoigt only supports matrices of dimension 1, 2 or 3");
        }
      },
      py::arg("matrix"), py::arg("isStrain") = true);

  /**
   * \brief Converts a vector given in Voigt notation to a matrix.
   *  \ingroup tensor
   * \param EVoigt Voigt notation vector.
   * \param isStrain Flag indicating whether the vector represents a strain (default is true).
   * \return Matrix corresponding to the vector in Voigt notation.
   *  \details
   * This function converts a vector given in Voigt notation to the corresponding matrix. The conversion depends on the
   * size The parameter `isStrain` is used to determine the conversion factor for off-diagonal components, which need to
   * be divided by 2 in the matrix representation if the quantity is a strain tensor.
   *
   * The function requires that the size of the Voigt notation vector is valid (1, 3, or 6).
   */
  utils.def(
      "fromVoigt",
      [](Eigen::MatrixXd vec, bool isStrain = true) {
        auto callFromVoigt = []<typename T>(const T& m, bool isStrain) {
          auto result = fromVoigt(m, isStrain);
          return Eigen::MatrixXd(result);
        };
        if (vec.rows() == 1) {
          Eigen::Vector<double, 1> fixedMat = vec;
          return callFromVoigt(fixedMat, isStrain);
        } else if (vec.rows() == 3) {
          Eigen::Vector<double, 3> fixedMat = vec;
          return callFromVoigt(fixedMat, isStrain);

        } else if (vec.rows() == 6) {
          Eigen::Vector<double, 6> fixedMat = vec;
          return callFromVoigt(fixedMat, isStrain);
        } else {
          DUNE_THROW(Dune::IOError, "fromVoigt only supports vectors of sizes 1, 3 or 6");
        }
      },
      py::arg("vector"), py::arg("isStrain") = true);
}