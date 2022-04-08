//
// Created by alex on 12/20/21.
//

#pragma once

#include "variableDefinitions.hh"

namespace Ikarus {
  enum class FEParameter { noParameter, loadfactor, time };
  enum class FESolutions { noSol, displacement, velocity, director, magnetizationAndVectorPotential };

  struct FEParameterValuePair {
    FEParameter type{FEParameter::noParameter};
    Ikarus::IVariable value;
  };

  struct FEParameterFactory {
    auto static createParameter(FEParameter&& feParameter, const int dimension) {
      switch (dimension) {
        case 1:
          return FEParameterValuePair(
              {feParameter, Ikarus::VariableFactory::createVariable(Ikarus::VariableTags::parameter1d)});
        case 2:
          return FEParameterValuePair(
              {feParameter, Ikarus::VariableFactory::createVariable(Ikarus::VariableTags::parameter2d)});
        case 3:
          return FEParameterValuePair(
              {feParameter, Ikarus::VariableFactory::createVariable(Ikarus::VariableTags::parameter3d)});
        default:
          assert(false && "You passed the wrong dimension size to FEParameterFactory::createParameter");
          __builtin_unreachable();
      }
    }
  };

}  // namespace Ikarus