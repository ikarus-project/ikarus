//
// Created by Alex on 19.05.2021.
//

#pragma once

#include "DefaultVariable.h"

#include "ikarus/Manifolds/RealTuple.h"
#include "ikarus/Manifolds/UnitVector.h"
#include <ikarus/Variables/GenericVariable.h>
#include <ikarus/Variables/VariableInterface.h>

#define IKARUS_REGISTER_VARIABLE(Manifold, ScalarType, Size, Tag, tagname) \
  using tagname = Variable<type<ScalarType, Size>, Tag> name;

namespace Ikarus::Variable {

  enum class VariablesTags : int {
    none,
    displacement1d,
    displacement2d,
    displacement3d,
    director2d,
    director3d,
    pressure
  };

  inline constexpr VariablesTags none = VariablesTags::none;
  inline constexpr VariablesTags displacement1d = VariablesTags::displacement1d;
  inline constexpr VariablesTags displacement2d = VariablesTags::displacement2d;
  inline constexpr VariablesTags displacement3d = VariablesTags::displacement3d;
  inline constexpr VariablesTags director2d = VariablesTags::director2d;
  inline constexpr VariablesTags director3d = VariablesTags::director3d;
  inline constexpr VariablesTags pressure = VariablesTags::pressure;

  using Manifold::RealTuple;
  using Manifold::UnitVector;

  using DISPLACEMENT1D = DefaultVariable<RealTuple<double, 1>, static_cast<int>(displacement1d)>;
  using DISPLACEMENT2D = DefaultVariable<RealTuple<double, 2>, static_cast<int>(displacement2d)>;
  using DISPLACEMENT3D = DefaultVariable<RealTuple<double, 3>, static_cast<int>(displacement3d)>;
  using DIRECTOR2D     = DefaultVariable<UnitVector<double, 2>, static_cast<int>(director2d)>;
  using DIRECTOR3D     = DefaultVariable<UnitVector<double, 3>, static_cast<int>(director3d)>;
  using PRESSURE       = DefaultVariable<UnitVector<double, 3>, static_cast<int>(pressure)>;

  class VariableFactory {
  public:
    static auto createVariable(VariablesTags tag) {
      switch (tag) {
        case displacement1d:
          return GenericVariable(DISPLACEMENT1D());
        case displacement2d:
          return GenericVariable(DISPLACEMENT2D());
        case displacement3d:
          return GenericVariable(DISPLACEMENT3D());
        case director2d:
          return GenericVariable(DIRECTOR2D());
        case director3d:
          return GenericVariable(DIRECTOR3D());
        case pressure:
          return GenericVariable(PRESSURE());
        default:
          throw std::logic_error("none variable not implemented. This should be an variable with empty size.");
      }
    }

  private:
    VariableFactory() = default;
  };


  //TODO: change this when gcc supports std::string,std::vector constexpr
  constexpr std::array<const char*,7> createVariableMap()
  {
    std::array<const char*,7> m{"none",
                                           "displacement1d",
                                           "displacement2d",
                                "displacement3d",
                                "director2d",
                                "director3d",
                                "pressure"    };

    return m;
  }

  inline constexpr std::array<const char*,7> variableNames = createVariableMap();


}  // namespace Ikarus::Variable
