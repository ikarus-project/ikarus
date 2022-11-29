// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "feRequirements.hh"

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
