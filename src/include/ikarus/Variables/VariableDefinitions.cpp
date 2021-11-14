//
// Created by Alex on 03.07.2021.
//
#include <ostream>

#include <ikarus/Variables/VariableDefinitions.h>
namespace Ikarus::Variable {
  std::ostream& operator<<(std::ostream& s, const VariableTags& varTag) {
    s << variableNames[static_cast<int>(varTag)];
    return s;
  }
}  // namespace Ikarus::Variable