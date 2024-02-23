// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <type_traits>

#include <ikarus/utils/tensorutils.hh>

/**
 *
 * \ingroup FEParameterTags
 * \brief A strongly typed enum class representing the type of the result request
 */

namespace Ikarus {
namespace Util {

  struct VectorizeWithVoigt
  {
    template <typename Derived>
    auto operator()(const Eigen::MatrixBase<Derived>& mat, const bool strainLike = false) const {
      return toVoigt(mat.derived(), strainLike);
    }
  };

  struct VectorizeGeneric
  {
    template <typename Derived>
    auto operator()(const Eigen::MatrixBase<Derived>& mat) const {
      return mat.derived().reshaped().eval();
    }
  };

  struct MatricizeWithVoigt
  {
    template <typename Derived>
    auto operator()(const Eigen::MatrixBase<Derived>& vec, const bool strainLike = false) const {
      return fromVoigt(vec.derived(), strainLike);
    }
  };

  struct MatricizeGeneric
  {
    template <typename Derived>
    auto operator()(const Eigen::MatrixBase<Derived>& vec) const {
      auto rc = static_cast<int>(std::floor(std::sqrt(vec.derived().size())));
      return vec.derived().reshaped(rc, rc).eval();
    }
  };

} // namespace Util

namespace ResultType {
#define REGISTER_RESULTTYPE_IMPL(structName, rowsExpr, ColsExpr, VectorizeStruct, MatricizeStruct) \
  template <typename ScalarType, int gridDim, int worldDim>                                        \
  struct structName                                                                                \
  {                                                                                                \
    friend auto toString(structName) { return #structName; }                                       \
                                                                                                   \
    using type       = Eigen::Matrix<ScalarType, rowsExpr, ColsExpr>;                              \
    using Vectorizer = VectorizeStruct;                                                            \
    using Matricizer = MatricizeStruct;                                                            \
  }

#define REGISTER_SYMMETRIC_RESULTTYPE(structName, rowsExpr, ColsExpr)                        \
  REGISTER_RESULTTYPE_IMPL(structName, rowsExpr, ColsExpr, Ikarus::Util::VectorizeWithVoigt, \
                           Ikarus::Util::MatricizeWithVoigt)

#define REGISTER_RESULTTYPE(structName, rowsExpr, ColsExpr)                                \
  REGISTER_RESULTTYPE_IMPL(structName, rowsExpr, ColsExpr, Ikarus::Util::VectorizeGeneric, \
                           Ikarus::Util::MatricizeGeneric)

  ///=============================
  /// Registered ResultTypes

  REGISTER_SYMMETRIC_RESULTTYPE(linearStress, worldDim, worldDim);
  REGISTER_SYMMETRIC_RESULTTYPE(PK2Stress, worldDim, worldDim);
  REGISTER_SYMMETRIC_RESULTTYPE(cauchyStress, worldDim, worldDim);

  REGISTER_RESULTTYPE(director, worldDim, 1);
  REGISTER_RESULTTYPE(magnetization, worldDim, 1);
  REGISTER_RESULTTYPE(gradientNormOfMagnetization, 1, 1);
  REGISTER_RESULTTYPE(vectorPotential, worldDim, 1);
  REGISTER_RESULTTYPE(divergenceOfVectorPotential, 1, 1);

  REGISTER_RESULTTYPE(BField, worldDim, 1);
  REGISTER_RESULTTYPE(HField, worldDim, 1);

  REGISTER_RESULTTYPE(customType, Eigen::Dynamic, Eigen::Dynamic);
} // namespace ResultType

template <template <typename, int, int> class RT1, template <typename, int, int> class RT2>
constexpr bool isSameResultType = std::is_same_v<RT1<double, 1, 1>, RT2<double, 1, 1>>;

/**
 * \brief Container that is used for FE Results. It gives access to the stored value, but can also be used to access
 * the result in Matrix and Vector form
 * \tparam RT A specified ResultType
 * \tparam inputIsVec boolean indicating wether the stored result is in its vector form, defaults to true
 */
template <typename RT, bool inputIsVec = true>
struct ResultTypeContainer : RT
{
private:
  using RegisteredType = typename RT::type;
  using Traits         = typename RegisteredType::CompileTimeTraits;
  using ScalarType     = typename RegisteredType::Scalar;
  static constexpr bool registeredTypeIsDynamic =
      Traits::RowsAtCompileTime == Eigen::Dynamic or Traits::ColsAtCompileTime == Eigen::Dynamic;
  static constexpr bool registeredTypeIsVector = Traits::IsVectorAtCompileTime;

public:
  using VecType    = std::conditional_t<registeredTypeIsDynamic, Eigen::MatrixX<ScalarType>,
                                     std::invoke_result_t<typename RT::Vectorizer, RegisteredType>>;
  using MatType    = std::conditional_t<registeredTypeIsDynamic, Eigen::MatrixX<ScalarType>,
                                     std::invoke_result_t<typename RT::Matricizer, VecType>>;
  using InputType  = std::conditional_t<inputIsVec, VecType, RegisteredType>;
  using ResultType = RT;

  const InputType& operator()() const { return value_; }

  /**
   * \brief Returns the stored value as Vector
   */
  auto asVec() const {
    if constexpr (inputIsVec)
      return value_;
    else
      return typename RT::Vectorizer{}(value_);
  }

  /**
   * \brief Returns the stored value as Matrix (if possible)
   */
  auto asMat() const {
    if constexpr (inputIsVec and not registeredTypeIsVector)
      return typename RT::Matricizer{}(value_);
    else
      return value_;
  }

  explicit ResultTypeContainer() = default;

  void emplace(const InputType&& value) { value_ = std::move(value); }

private:
  InputType value_{};
};

template <template <typename, int, int> class RT>
auto toString() {
  return toString(RT<double, 1, 1>{});
}

} // namespace Ikarus