//
// Created by Alex on 16.06.2021.
//

#include "finiteElementFunctionConcepts.hh"

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
      case ResultType::divergenceOfVectorPotential:
        return "divergenceOfVectorPotential";
        break;
      case ResultType::cauchyStress:
        return "cauchyStress";
        break;
      case ResultType::director:
        return "director";
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
