//
// Created by Alex on 19.05.2021.
//

#pragma once

#include "defaultVariable.hh"

#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/manifolds/unitVector.hh>
#include <ikarus/variables/interfaceVariable.hh>
#include <ikarus/variables/variablePolicies.hh>

namespace Ikarus {

  // These can be used as variables of some fe solver, i.e. "Degrees of Freedom" or just as data
  enum class VariableTags : int {
    displacement1d,
    displacement2d,
    displacement3d,
    director2d,
    director3d,
    pressure,
    velocity1d,
    velocity2d,
    edgeLength,
    parameter1d,
    parameter2d,
    parameter3d
  };

  inline constexpr VariableTags AllVariableTags[]
      = {VariableTags::displacement1d, VariableTags::displacement2d, VariableTags::displacement3d,
         VariableTags::director2d,     VariableTags::director3d,     VariableTags::pressure,
         VariableTags::velocity1d,     VariableTags::velocity2d,     VariableTags::edgeLength,
         VariableTags::parameter1d,    VariableTags::parameter2d,    VariableTags::parameter3d};

  using Manifold::RealTuple;
  using Manifold::UnitVector;

  using DISPLACEMENT1D = DefaultVariable<RealTuple<double, 1>, static_cast<int>(VariableTags::displacement1d)>;
  using DISPLACEMENT2D = DefaultVariable<RealTuple<double, 2>, static_cast<int>(VariableTags::displacement2d)>;
  using DISPLACEMENT3D = DefaultVariable<RealTuple<double, 3>, static_cast<int>(VariableTags::displacement3d)>;
  using DIRECTOR2D     = DefaultVariable<UnitVector<double, 2>, static_cast<int>(VariableTags::director2d)>;
  using DIRECTOR3D     = DefaultVariable<UnitVector<double, 3>, static_cast<int>(VariableTags::director3d)>;
  using PRESSURE       = DefaultVariable<RealTuple<double, 1>, static_cast<int>(VariableTags::pressure)>;
  using VELOCITY1D     = DefaultVariable<RealTuple<double, 1>, static_cast<int>(VariableTags::velocity1d)>;
  using VELOCITY2D     = DefaultVariable<RealTuple<double, 2>, static_cast<int>(VariableTags::velocity2d)>;
  using EDGELENGTH     = DefaultVariable<RealTuple<double, 1>, static_cast<int>(VariableTags::edgeLength)>;
  using PARAMETER1D    = DefaultVariable<RealTuple<double, 1>, static_cast<int>(VariableTags::parameter1d)>;
  using PARAMETER2D    = DefaultVariable<RealTuple<double, 2>, static_cast<int>(VariableTags::parameter2d)>;
  using PARAMETER3D    = DefaultVariable<RealTuple<double, 3>, static_cast<int>(VariableTags::parameter3d)>;

  class VariableFactory {
  public:
    static auto createVariable(VariableTags tag) {
      switch (tag) {
        case VariableTags::displacement1d:
          return IVariable(DISPLACEMENT1D());
        case VariableTags::displacement2d:
          return IVariable(DISPLACEMENT2D());
        case VariableTags::displacement3d:
          return IVariable(DISPLACEMENT3D());
        case VariableTags::director2d:
          return IVariable(DIRECTOR2D());
        case VariableTags::director3d:
          return IVariable(DIRECTOR3D());
        case VariableTags::pressure:
          return IVariable(PRESSURE());
        case VariableTags::velocity1d:
          return IVariable(VELOCITY1D());
        case VariableTags::velocity2d:
          return IVariable(VELOCITY2D());
        case VariableTags::edgeLength:
          return IVariable(EDGELENGTH());
        case VariableTags::parameter1d:
          return IVariable(PARAMETER1D());
        case VariableTags::parameter2d:
          return IVariable(PARAMETER2D());
        case VariableTags::parameter3d:
          return IVariable(PARAMETER3D());
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

  std::ostream& operator<<(std::ostream& s, const VariableTags& varTag);

}  // namespace Ikarus