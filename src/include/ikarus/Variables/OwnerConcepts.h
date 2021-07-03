//
// Created by Alex on 26.05.2021.
//

#pragma once

#include <ikarus/Variables/VariablePolicies.h>

namespace Ikarus::Concepts {
  template <typename VariableOwnerType>
  concept VariableOwner = requires(VariableOwnerType varOwner, Variable var) {
    typename Node::ctype;
    varOwner.addVariable(var);
    varOwner.removeVariable(var);
    varOwner.getVariable(var);
    varOwner.setVariable(var);
  };

  template <typename DOFOwnerType>
  concept DOFOwner = requires(DOFOwnerType dofOwner, Variable var) {
    typename DOFOwner::ctype;
    dofOwner.addDOF(var);
    dofOwner.removeDOF(var);
    dofOwner.getDOF(var);
    dofOwner.setDOF(var);
  };

};  // namespace Ikarus::Concepts