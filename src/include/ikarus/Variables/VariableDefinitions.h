//
// Created by Alex on 19.05.2021.
//

#pragma once

#include "DefaultVariable.h"

#include "ikarus/Manifolds/RealTuple.h"
#include "ikarus/Manifolds/UnitVector.h"
#include <ikarus/Variables/VariableInterface.h>

#define IKARUS_DEFINE_VECTOR_VARIABLE(Manifold, ScalarType, Size, Tag, tagname) \
  using tagname = Variable<type<ScalarType, Size>, Tag> name;

namespace Ikarus::Variable {

  template <size_t IntTag> struct VariableTag {
    /** \brief Tag to distinguish between variable definitions */
    static constexpr size_t tagValue = IntTag;
  };
  template <size_t tagval1, size_t tagval2> bool operator<(VariableTag<tagval1>, VariableTag<tagval2>) {
    return tagval1 < tagval2;
  }

  template <size_t tagval1, size_t tagval2> bool operator==(VariableTag<tagval1>, VariableTag<tagval2>) {
    return tagval1 == tagval2;
  }

  template <size_t tagval1, size_t tagval2> bool operator>(VariableTag<tagval1>, VariableTag<tagval2>) {
    return tagval1 > tagval2;
  }

  template <size_t tagval1, size_t tagval2> bool operator!=(VariableTag<tagval1>, VariableTag<tagval2>) {
    return tagval1 != tagval2;
  }

  template <size_t tagval1, size_t tagval2> bool operator>=(VariableTag<tagval1>, VariableTag<tagval2>) {
    return tagval1 >= tagval2;
  }

  template <size_t tagval1, size_t tagval2> bool operator<=(VariableTag<tagval1>, VariableTag<tagval2>) {
    return tagval1 <= tagval2;
  }

  using Manifold::RealTuple;
  using Manifold::UnitVector;
  using DISPLACEMENT1D = DefaultVariable<RealTuple<double, 1>, VariableTag<0>>;
  using DISPLACEMENT2D = DefaultVariable<RealTuple<double, 2>, VariableTag<1>>;
  using DISPLACEMENT3D = DefaultVariable<RealTuple<double, 3>, VariableTag<2>>;
  using DIRECTOR2D = DefaultVariable<UnitVector<double, 3>, VariableTag<3>>;
  using DIRECTOR3D = DefaultVariable<UnitVector<double, 3>, VariableTag<4>>;
}  // namespace Ikarus::Variable
