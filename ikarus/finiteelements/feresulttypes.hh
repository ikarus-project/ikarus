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
namespace Impl {

  template<bool strainlike=false>
  struct VectorizeWithVoigt
  {
    template <typename Derived>
    static auto transform(const Eigen::DenseBase<Derived>& mat)  {
      return toVoigt(mat.derived(), strainlike);
    }
  };

  struct VectorizeGeneric
  {
    template <typename Derived>
    static auto transform(const Eigen::DenseBase<Derived>& mat) {
      return mat.derived().reshaped().eval();
    }
  };

  template<bool strainlike=false>
  struct MatricizeWithVoigt
  {
    template <typename Derived, int RowsAtCompileTime, int ColsAtCompileTime>
    static auto transform(const Eigen::DenseBase<Derived>& vec,int rows= RowsAtCompileTime, int cols= ColsAtCompileTime) {
      assert(rows==RowsAtCompileTime && cols== ColsAtCompileTime && "Only the fixed size values work for voigt matrices and vectors");
      static_assert(RowsAtCompileTime!= Eigen::Dynamic and ColsAtCompileTime!= Eigen::Dynamic,"Voigt notation only available for fixed size vectors and matrices");
      return fromVoigt(vec.derived(), strainlike);
    }
  };

  struct MatricizeGeneric
  {
    template <typename Derived, int RowsAtCompileTime, int ColsAtCompileTime>
    static auto transform(const Eigen::DenseBase<Derived>& vec,int rows= RowsAtCompileTime, int cols= ColsAtCompileTime) {
      return vec.derived().reshaped(Eigen::fix<RowsAtCompileTime>(rows), Eigen::fix<ColsAtCompileTime>(cols)).eval();
    }


  };

} // namespace Impl

namespace ResultType {
#define REGISTER_RESULTTYPE_IMPL(structName, rowsExpr, ColsExpr,MaxRowsExpr,MaxColsExpr, VectorizeStruct, MatricizeStruct) \
  template <typename ScalarType, int gridDim, int worldDim>                                        \
  struct structName                                                                                \
  {                                                                                                \
    friend auto toString(structName) { return #structName; }                                       \
                                                                                                   \
    using type       = Eigen::Matrix<ScalarType, rowsExpr, ColsExpr,0,MaxRowsExpr,MaxColsExpr>;                              \
    using Vectorizer = VectorizeStruct;                                                            \
    using Matricizer = MatricizeStruct;                                                            \
  }

#define REGISTER_SYMMETRIC_RESULTTYPE(structName, rowsExpr, ColsExpr,MaxRowsExpr,MaxColsExpr)                        \
  REGISTER_RESULTTYPE_IMPL(structName, rowsExpr, ColsExpr, Ikarus::Util::VectorizeWithVoigt, \
                           Ikarus::Util::MatricizeWithVoigt)

  #define REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(structName, rowsExpr, ColsExpr,strainlike)                        \
  REGISTER_RESULTTYPE_IMPL(structName, rowsExpr, ColsExpr, rowsExpr, ColsExpr, Ikarus::Util::VectorizeWithVoigt<strainlike>, \
  Ikarus::Util::MatricizeWithVoigt<strainlike>)

#define REGISTER_RESULTTYPE(structName, rowsExpr, ColsExpr,MaxRowsExpr,MaxColsExpr)                                \
  REGISTER_RESULTTYPE_IMPL(structName, rowsExpr, ColsExpr, Ikarus::Util::VectorizeGeneric, \
                           Ikarus::Util::MatricizeGeneric)

  #define REGISTER_RESERVED_RESULTTYPE(structName, rowsExpr, ColsExpr,MaxRowsExpr,MaxColsExpr)                                \
REGISTER_RESULTTYPE_IMPL(structName, rowsExpr, ColsExpr, MaxRowsExpr,MaxColsExpr, Ikarus::Util::VectorizeGeneric, \
Ikarus::Util::MatricizeGeneric)

  #define REGISTER_SIMPLE_RESULTTYPE(structName, rowsExpr, ColsExpr)                                \
REGISTER_RESERVED_RESULTTYPE(structName, rowsExpr, ColsExpr, rowsExpr, ColsExpr)

  ///=============================
  /// Registered ResultTypes

  REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(linearStress, worldDim, worldDim,false);
  REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(PK2Stress, worldDim, worldDim,false);
  REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(cauchyStress, worldDim, worldDim,false);

  REGISTER_SIMPLE_RESULTTYPE(director, worldDim, 1);
  REGISTER_SIMPLE_RESULTTYPE(magnetization, worldDim, 1);
  REGISTER_SIMPLE_RESULTTYPE(gradientNormOfMagnetization, 1, 1);
  REGISTER_SIMPLE_RESULTTYPE(vectorPotential, worldDim, 1);
  REGISTER_SIMPLE_RESULTTYPE(divergenceOfVectorPotential, 1, 1);

  REGISTER_SIMPLE_RESULTTYPE(BField, worldDim, 1);
  REGISTER_SIMPLE_RESULTTYPE(HField, worldDim, 1);

  REGISTER_SIMPLE_RESULTTYPE(customType, Eigen::Dynamic, Eigen::Dynamic);
} // namespace ResultType

template <template <typename, int, int> class RT1, template <typename, int, int> class RT2>
constexpr bool isSameResultType = std::is_same_v<RT1<double, 1, 1>, RT2<double, 1, 1>>;

/**
 * \brief Container that is used for FE Results. It gives access to the stored value, but can also be used to access
 * the result in Matrix and Vector form
 * \tparam RT A specified ResultType
 * \tparam storedValueIsVector boolean indicating wether the stored result is in its vector form, defaults to true
 */
template <typename RT, bool storedValueIsVector=true>
struct ResultTypeContainer : RT
{
private:
  using ResultTypeValueType = typename RT::type;
  static constexpr Eigen::Index rowsAtCompileTime = ResultTypeValueType::RowsAtCompileTime;
  static constexpr Eigen::Index colsAtCompileTime = ResultTypeValueType::ColsAtCompileTime;

public:
  using VecType    = std::invoke_result_t<decltype(&RT::Vectorizer::template transform<ResultTypeValueType>), const ResultTypeValueType&>;
  using MatType    =  std::invoke_result_t<decltype(&RT::Matricizer::template transform<VecType,rowsAtCompileTime,colsAtCompileTime>),const VecType&,int,int>;
  using StoredType = std::conditional_t<storedValueIsVector,VecType,MatType>;
  using ResultType = RT;

  /**
   * \brief Returns the stored value as Vector
   */
  auto asVec() const {
    if constexpr (storedValueIsVector)
      return value_;
    else
      return RT::Vectorizer::transform(value_);
  }

  /**
   * \brief Returns the stored value as Matrix (if possible)
   */

  auto asMat(Eigen::Index rows =rowsAtCompileTime, Eigen::Index cols =colsAtCompileTime) const {
    if constexpr (storedValueIsVector) {
      if constexpr (rowsAtCompileTime==Eigen::Dynamic )
      assert( rows!=rowsAtCompileTime  &&
        "For dynamic size result types you have to pass rows by hand, since it is not clear how this should be reshaped" );
      if constexpr (colsAtCompileTime==Eigen::Dynamic )
        assert( cols!=colsAtCompileTime &&
          "For dynamic size result types you have to pass cols by hand, since it is not clear how this should be reshaped" );
      return RT:: Matricizer::template transform<VecType,rowsAtCompileTime,colsAtCompileTime>(value_,rows,cols);
    }
    else
      return value_;
  }
  explicit ResultTypeContainer() = default;

  explicit ResultTypeContainer(StoredType&& value) {this->value_=std::move(value); }
  explicit ResultTypeContainer(const StoredType& value) {this->value_=value;}
  ResultTypeContainer& operator=(const StoredType& value){this->value_=value;  return *this;}
  ResultTypeContainer& operator=( StoredType&& value){this->value_=std::move(value);  return *this;}

private:
  StoredType value_{};
};

template <template <typename, int, int> class RT>
auto toString() {
  return toString(RT<double, 1, 1>{});
}

} // namespace Ikarus
