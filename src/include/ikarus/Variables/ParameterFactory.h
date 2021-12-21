//
// Created by alex on 12/20/21.
//

#pragma once

#include "VariableDefinitions.h"

namespace Ikarus {
  enum class FEParameter { noParameter, loadfactor };

  struct FEParameterValuePair {
    FEParameter type{FEParameter::noParameter};
    Ikarus::Variable::IVariable value;
  };

  struct FEParameterFactory {
    auto static createParameter(FEParameter&& feParameter, const int dimension) {
      switch (dimension) {
        case 1:
          return FEParameterValuePair({feParameter, Ikarus::Variable::VariableFactory::createVariable(
                                                        Ikarus::Variable::VariableTags::parameter1d)});
        case 2:
          return FEParameterValuePair({feParameter, Ikarus::Variable::VariableFactory::createVariable(
                                                        Ikarus::Variable::VariableTags::parameter2d)});
        case 3:
          return FEParameterValuePair({feParameter, Ikarus::Variable::VariableFactory::createVariable(
                                                        Ikarus::Variable::VariableTags::parameter3d)});
        default:
          assert(false && "You passed the wrong dimension size to FEParameterFactory::createParameter");
          __builtin_unreachable();
      }
    }
  };

}  // namespace Ikarus