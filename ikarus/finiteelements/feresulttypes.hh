// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/utils/math.hh>
/**
 *
 * \ingroup FEParameterTags
 * \brief A strongly typed enum class representing the type of the result request
 */

namespace Ikarus::ResultType {
#define REGISTER_RT(structName) \
  friend auto toString(structName) { return #structName; }

namespace Util {
  template <int dim>
  constexpr int voigtSize() {
    return dim * (dim + 1) / 2;
  }

  template <int dim>
  constexpr int matrixSize() {
    return (-1 + ct_sqrt(1 + 8 * voigtSize<dim>())) / 2;
  }

} // namespace Util

struct noType
{
  REGISTER_RT(noType);
};

struct linearStress
{
  REGISTER_RT(linearStress);

  using voigtApplicable = std::true_type;

  template <int dim>
  using matrixType = Eigen::Matrix<double, Util::matrixSize<dim>(), Util::matrixSize<dim>()>;

  template <int dim>
  using voigtType = Eigen::Vector<double, Util::voigtSize<dim>()>;

  template <int gridDim, int worldDim, bool voigt>
  using type = std::conditional_t<voigt, voigtType<gridDim>, matrixType<gridDim>>;
};

struct PK2Stress
{
  REGISTER_RT(PK2Stress);

  using voigtApplicable = std::true_type;

  template <int dim>
  using matrixType = Eigen::Matrix<double, Util::matrixSize<dim>(), Util::matrixSize<dim>()>;

  template <int dim>
  using voigtType = Eigen::Vector<double, Util::voigtSize<dim>()>;

  template <int gridDim, int worldDim, bool voigt>
  using type = std::conditional_t<voigt, voigtType<gridDim>, matrixType<gridDim>>;
};

struct cauchyStress
{
  REGISTER_RT(cauchyStress);

  using voigtApplicable = std::true_type;

  template <int dim>
  using matrixType = Eigen::Matrix<double, Util::matrixSize<dim>(), Util::matrixSize<dim>()>;

  template <int dim>
  using voigtType = Eigen::Vector<double, Util::voigtSize<dim>()>;

  template <int gridDim, int worldDim, bool voigt>
  using type = std::conditional_t<voigt, voigtType<gridDim>, matrixType<gridDim>>;
};

struct director
{
  REGISTER_RT(director);

  template <int gridDim, int worldDim>
  using type = Eigen::Vector<double, worldDim>;
};

struct magnetization
{
  REGISTER_RT(magnetization);
};

struct gradientNormOfMagnetization
{
  REGISTER_RT(gradientNormOfMagnetization);
};

struct vectorPotential
{
  REGISTER_RT(vectorPotential);
};

struct divergenceOfVectorPotential
{
  REGISTER_RT(divergenceOfVectorPotential);
};

struct BField
{
  REGISTER_RT(BField);
};

struct HField
{
  REGISTER_RT(HField);
};

struct customType
{
  REGISTER_RT(customType);

  template <int gridDim, int worldDim>
  using type = Eigen::MatrixXd;
};

template <typename RT>
concept ResultTypeConcept = requires(RT t) {
  { toString(t) } -> std::convertible_to<std::string>;
};

template <typename RT>
concept HasVoigt = std::is_same_v<typename RT::voigtApplicable, std::true_type> && ResultTypeConcept<RT>;

} // namespace Ikarus::ResultType

namespace Ikarus {

template <typename RT, int gridDim, int worldDim, bool>
struct getResultType
{
  using type = typename RT::template type<gridDim, worldDim>;
};

template <ResultType::HasVoigt RT, int gridDim, int worldDim, bool voigt>
struct getResultType<RT, gridDim, worldDim, voigt>
{
  using type = typename RT::template type<gridDim, worldDim, voigt>;
};

template <typename RT, int gridDim, int worldDim = gridDim, bool voigt = true>
using resultType_t = typename getResultType<RT, gridDim, worldDim, voigt>::type;

} // namespace Ikarus