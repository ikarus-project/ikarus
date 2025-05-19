// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#if ENABLE_MUESLI

  #include <muesli/muesli.h>

  #include <ikarus/finiteelements/physicshelper.hh>
  #include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Concepts {
/**
 * \concept MuesliMaterialImplementation
 * \brief Concept to check weather a type is a Mueli material implementation
 *
 * \tparam MAT the material type
 */
template <typename MAT>
concept MuesliMaterialImplementation = std::is_base_of_v<muesli::material, MAT>;
} // namespace Ikarus::Concepts

namespace Ikarus::Materials {

/**
 * \brief converts Ikarus material parameters to Muesli material properties.
 *
 * \tparam MPT the type of the Ikarus material parameters.
 * \param mpt the Ikarus material parameters.
 * \return MaterialProperties muesli material properties with the Lambda and mu set
 */
template <Concepts::MPTuple MPT>
inline muesli::materialProperties propertiesFromIkarusMaterialParameters(const MPT& mpt) {
  auto converter = convertLameConstants(mpt);

  auto mpm = muesli::materialProperties{};
  mpm.insert({"lambda", converter.toLamesFirstParameter()});
  mpm.insert({"mu", converter.toShearModulus()});

  return mpm;
}

/**
 * \brief adds a specific tag to the muesli materialproperties, thic can be used for example to add the `regularized`
 * tag for Neo-Hooke or the `compressible` tag for Yeoh and Arruda-Boyce.
 */
inline void addTag(muesli::materialProperties& mpm, const std::string& tagName, double tagValue = 0) {
  mpm.insert({tagName, tagValue});
}

/**
 * \brief Converts the entries of a Eigen::Matrix to a provided muesli::istensor (symmetric 2nd order tensor).
 *
 * \param C the Eigen::Matrix that is to be converted.
 */
inline istensor toistensor(const Eigen::Matrix<double, 3, 3>& C) {
  return istensor(C(0, 0), C(1, 1), C(2, 2), C(1, 2), C(2, 0), C(0, 1));
}

/**
 * \brief Converts the entries of a Eigen::Matrix to a provided muesli::itensor (2nd order tensor).
 *
 * \param C the Eigen::Matrix that is to be converted.
 */
inline itensor toitensor(const Eigen::Matrix<double, 3, 3>& C) {
  return itensor(C(0, 0), C(0, 1), C(0, 2), C(1, 0), C(1, 1), C(1, 2), C(2, 0), C(2, 1), C(2, 2));
}

/**
 * \brief Converts a provided muesli::istensor (symmetric 2nd order tensor) to a Eigen::Matrix.
 *
 * \param it provided istensor.
 * \return Eigen::Matrix<ScalarType, dim, dim> converted matrix.
 */
inline auto toEigenMatrix(const istensor& it) {
  auto S = Eigen::Matrix<double, 3, 3>{};
  for (auto i : Dune::range(3))
    for (auto j : Dune::range(3))
      S(i, j) = it(i, j);
  return S;
}

/**
 * \brief Converts a provided muesli::istensor4 (symmetric 4th order tensor) to a Eigen::TensorFixedSize.
 *
 * \param it provided itensor.
 * \return Eigen::TensorFixedSize<ScalarType, Eigen::Sizes<dim, dim, dim, dim>> converted tensor.
 */
inline auto toEigenTensor(const itensor4& it) {
  Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>> moduli{};
  for (auto i : Dune::range(3))
    for (auto j : Dune::range(3))
      for (auto k : Dune::range(3))
        for (auto l : Dune::range(3))
          moduli(i, j, k, l) = it(i, j, k, l);

  return moduli;
}

/**
 * \brief Generates the name of the Muesli material from the class name.
 *
 * \tparam MAT type of the Muesli material implementation.
 * \return constexpr std::string the name of the material.
 */
template <Concepts::MuesliMaterialImplementation MAT>
constexpr std::string materialName() {
  std::string matName = Dune::className<MAT>();

  // Find and remove the "muesli::" namespace
  std::string prefix = "muesli::";
  matName            = matName.substr(prefix.size());

  // Capitalize the first letter
  if (!matName.empty() && std::islower(matName[0])) {
    matName[0] = std::toupper(matName[0]);
  }
  return matName;
}

} // namespace Ikarus::Materials

#else
  #error Muesli materials depends on the Muesli library, which is not included
#endif