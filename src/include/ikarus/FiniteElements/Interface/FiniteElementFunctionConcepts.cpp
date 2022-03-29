//
// Created by Alex on 16.06.2021.
//

#include "FiniteElementFunctionConcepts.h"

namespace Ikarus {

  std::string getResultType(const ResultType& res) {
    switch (res) {
      case ResultType::noType:
        return "noType";
      case ResultType::magnetization:
        return "magnetization";
      case ResultType::vectorPotential:
        return "vectorPotential";
      case ResultType::gradientNormOfMagnetization:
        return "gradientNormOfMagnetization";
        break;
      case ResultType::BField:
        return "BField";
        break;
      case ResultType::HField:
        return "HField";
        break;
    }
    __builtin_unreachable();
  }

}  // namespace Ikarus
