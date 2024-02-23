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
    auto operator()(const Eigen::MatrixBase<Derived>& mat) const {
      return toVoigt(mat.derived(), false);
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
    auto operator()(const Eigen::MatrixBase<Derived>& vec) const {
      return fromVoigt(vec.derived(), false);
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

  REGISTER_RESULTTYPE(BField, 1, 1);
  REGISTER_RESULTTYPE(HField, 1, 1);

  REGISTER_RESULTTYPE(customType, Eigen::Dynamic, Eigen::Dynamic);
} // namespace ResultType

template <template <typename, int, int> class RT1, template <typename, int, int> class RT2>
constexpr bool isSameResultType = std::is_same_v<RT1<double, 1, 1>, RT2<double, 1, 1>>;

template <typename RT, bool inputIsVec>
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
  const InputType& operator()() { return value_; }

  auto asVec() const {
    if constexpr (inputIsVec)
      return value_;
    else
      return vectorizer_(value_);
  }
  auto asMat() const {
    if constexpr (inputIsVec and not registeredTypeIsVector)
      return matricizer_(value_);
    else
      return value_;
  }

  void emplace(const InputType& t) { value_ = std::move(t); }

  explicit ResultTypeContainer() = default;

private:
  InputType value_{};
  typename RT::Vectorizer vectorizer_{};
  typename RT::Matricizer matricizer_{};
};

template <template <typename, int, int> class RT>
auto toString() {
  return toString(RT<double, 1, 1>{});
}

} // namespace Ikarus