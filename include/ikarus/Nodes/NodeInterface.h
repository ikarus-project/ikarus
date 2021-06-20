//
// Created by Alex on 19.05.2021.
//

#ifndef IKARUS_NODEINTERFACE_H
#define IKARUS_NODEINTERFACE_H

#include <ikarus/Variables/OwnerConcepts.h>

namespace Ikarus::Concepts {
  template <typename NodeType>
  concept Node
      = Ikarus::VariableOwner && Ikarus::VariableOwner && requires(NodeType node, Variable var) {
    typename Node::ctype;
  };

};  // namespace Ikarus::Concepts

#endif  // IKARUS_NODEINTERFACE_H
