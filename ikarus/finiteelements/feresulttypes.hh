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

namespace Impl {
  template <int dim>
  constexpr int voigtSize() {
    return dim * (dim + 1) / 2;
  }

  template <int dim>
  constexpr int matrixSize() {
    return (-1 + ct_sqrt(1 + 8 * voigtSize<dim>())) / 2;
  }

} // namespace Impl

struct noType;
struct magnetization;
struct gradientNormOfMagnetization;
struct vectorPotential;
struct divergenceOfVectorPotential;
struct BField;
struct HField;
struct cauchyStress;
struct PK2Stress;
struct linearStress;
struct director;
struct customType;

struct noType
{
  REGISTER_RT(noType);
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
struct cauchyStress
{
  REGISTER_RT(cauchyStress);

  static constexpr bool voigtApplicable = true;

  template <int gridDim>
  using matrixType = Eigen::Matrix<double, Impl::matrixSize<gridDim>(), Impl::matrixSize<gridDim>()>;

  template <int gridDim>
  using voigtType = Eigen::Vector<double, Impl::voigtSize<gridDim>()>;

  template <int gridDim, bool voigt>
  using type = std::conditional_t<voigt, voigtType<gridDim>, matrixType<gridDim>>;
};
struct PK2Stress
{
  REGISTER_RT(PK2Stress);

  static constexpr bool voigtApplicable = true;

  template <int gridDim>
  using matrixType = Eigen::Matrix<double, Impl::matrixSize<gridDim>(), Impl::matrixSize<gridDim>()>;

  template <int gridDim>
  using voigtType = Eigen::Vector<double, Impl::voigtSize<gridDim>()>;

  template <int gridDim, bool voigt>
  using type = std::conditional_t<voigt, voigtType<gridDim>, matrixType<gridDim>>;
};
struct linearStress
{
  REGISTER_RT(linearStress);

  static constexpr bool voigtApplicable = true;

  template <int gridDim>
  using matrixType = Eigen::Matrix<double, Impl::matrixSize<gridDim>(), Impl::matrixSize<gridDim>()>;

  template <int gridDim>
  using voigtType = Eigen::Vector<double, Impl::voigtSize<gridDim>()>;

  template <int gridDim, bool voigt>
  using type = std::conditional_t<voigt, voigtType<gridDim>, matrixType<gridDim>>;
};
struct director
{
  REGISTER_RT(director);

  static constexpr bool voigtApplicable = false;

  template <int gridDim>
  using type = Eigen::Vector<double, gridDim>;
};

struct customType
{
  REGISTER_RT(customType);

  template <int gridDim>
  using type = Eigen::MatrixXd;
};

} // namespace Ikarus::ResultType

namespace Ikarus {

template <typename RT>
concept HasVoigt = std::same_as<RT, ResultType::linearStress> || std::same_as<RT, ResultType::cauchyStress> ||
                   std::same_as<RT, ResultType::PK2Stress>;

template <typename RT, int gridDim, bool voigt>
struct getResultType
{
  using type = typename RT::template type<gridDim>;
};

template <HasVoigt RT, int gridDim, bool voigt>
struct getResultType<RT, gridDim, voigt>
{
  using type = typename RT::template type<gridDim, voigt>;
};

template <typename RT, int gridDim, bool voigt = true>
using resultType_t = typename getResultType<RT, gridDim, voigt>::type;

} // namespace Ikarus