//
// Created by Alex on 19.05.2021.
//

#ifndef IKARUS_DOFHANDLERINTERFACE_H
#define IKARUS_DOFHANDLERINTERFACE_H

#include "OccupationPattern.h"
#include "ikarus/Variables/VariableInterface.h"

namespace Ikarus {
  template <typename DofHandlerType>
  concept DofHandlerConcept
      = requires(DofHandlerType dh, GenericVariable var, GenericVariableOwner d) {
    dh.generateOccupationPattern();
    { dh.getOccupationPattern() } -> std::same_as<OccupationPattern>;
    { dh.getSolution(d, var) } -> std::same_as<OccupationPattern>;
  };
}  // namespace Ikarus
#endif  // IKARUS_DOFHANDLERINTERFACE_H
