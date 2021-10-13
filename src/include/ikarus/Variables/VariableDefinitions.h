//
// Created by Alex on 19.05.2021.
//

#pragma once

#include "DefaultVariable.h"

#include "ikarus/Manifolds/RealTuple.h"
#include "ikarus/Manifolds/UnitVector.h"
#include <ikarus/Variables/InterfaceVariable.h>
#include <ikarus/Variables/VariablePolicies.h>

namespace Ikarus::Variable {

  // These can be used as variables of some fe solver, i.e. "Degrees of Freedom" or just as data
  enum class VariablesTags : int {
    displacement1d,
    displacement2d,
    displacement3d,
    director2d,
    director3d,
    pressure,
    velocity1d,
    velocity2d,
    edgeLength
  };

  inline constexpr VariablesTags displacement1d = VariablesTags::displacement1d;
  inline constexpr VariablesTags displacement2d = VariablesTags::displacement2d;
  inline constexpr VariablesTags displacement3d = VariablesTags::displacement3d;
  inline constexpr VariablesTags director2d     = VariablesTags::director2d;
  inline constexpr VariablesTags director3d     = VariablesTags::director3d;
  inline constexpr VariablesTags pressure       = VariablesTags::pressure;
  inline constexpr VariablesTags velocity1d     = VariablesTags::velocity1d;
  inline constexpr VariablesTags velocity2d     = VariablesTags::velocity2d;
  inline constexpr VariablesTags edgeLength     = VariablesTags::edgeLength;

  inline constexpr VariablesTags AllVariableTags[]
      = {displacement1d, displacement2d, displacement3d, director2d, director3d,
         pressure,       velocity1d,     velocity2d,     edgeLength};

  using Manifold::RealTuple;
  using Manifold::UnitVector;

  using DISPLACEMENT1D = DefaultVariable<RealTuple<double, 1>, static_cast<int>(displacement1d)>;
  using DISPLACEMENT2D = DefaultVariable<RealTuple<double, 2>, static_cast<int>(displacement2d)>;
  using DISPLACEMENT3D = DefaultVariable<RealTuple<double, 3>, static_cast<int>(displacement3d)>;
  using DIRECTOR2D     = DefaultVariable<UnitVector<double, 2>, static_cast<int>(director2d)>;
  using DIRECTOR3D     = DefaultVariable<UnitVector<double, 3>, static_cast<int>(director3d)>;
  using PRESSURE       = DefaultVariable<RealTuple<double, 1>, static_cast<int>(pressure)>;
  using VELOCITY1D     = DefaultVariable<RealTuple<double, 1>, static_cast<int>(velocity1d)>;
  using VELOCITY2D     = DefaultVariable<RealTuple<double, 2>, static_cast<int>(velocity2d)>;
  using EDGELENGTH     = DefaultVariable<RealTuple<double, 1>, static_cast<int>(edgeLength)>;

  class VariableFactory {
  public:
    static auto createVariable(VariablesTags tag) {
      switch (tag) {
        case displacement1d:
          return IVariable(DISPLACEMENT1D());
        case displacement2d:
          return IVariable(DISPLACEMENT2D());
        case displacement3d:
          return IVariable(DISPLACEMENT3D());
        case director2d:
          return IVariable(DIRECTOR2D());
        case director3d:
          return IVariable(DIRECTOR3D());
        case pressure:
          return IVariable(PRESSURE());
        case velocity1d:
          return IVariable(VELOCITY1D());
        case velocity2d:
          return IVariable(VELOCITY2D());
        case edgeLength:
          return IVariable(EDGELENGTH());
        default:
          throw std::logic_error("Variable not implemented.");
      }
    }

  private:
    VariableFactory() = default;
  };

  // TODO: change this when gcc supports std::string,std::vector constexpr
  constexpr std::array<const char*, 9> createVariableMap() {
    std::array<const char*, 9> m{"displacement1d", "displacement2d", "displacement3d", "director2d", "director3d",
                                 "pressure",       "velocity1d",     "velocity2d",     "edgeLength"};

    return m;
  }

  inline constexpr std::array<const char*, 9> variableNames = createVariableMap();

  std::ostream& operator<<(std::ostream& s, const VariablesTags& varTag);

}  // namespace Ikarus::Variable