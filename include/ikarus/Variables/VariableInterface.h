//
// Created by Alex on 19.05.2021.
//

#ifndef IKARUS_VARIABLE_INTERFACE_H
#define IKARUS_VARIABLE_INTERFACE_H

#include <concepts>
#include <numeric>

#include "dune/common/classname.hh"

namespace Ikarus::Concepts {

  /**
   * \brief Concept of a variable
   *
   *
   */
  template <typename VariableType>
  concept Variable = requires(VariableType var, VariableType var2,
                              typename VariableType::CorrectionType correction,
                              typename VariableType::CoordinateType value) {
    typename VariableType::ctype;
    VariableType::valueSize;
    VariableType::correctionSize;
    typename VariableType::CoordinateType;
    typename VariableType::CorrectionType;
    { var.getValue() } -> std::convertible_to<typename VariableType::CoordinateType>;
    { var.setValue(value) } -> std::same_as<void>;
    { var.update(correction) } -> std::same_as<void>;
    { var == var2 } -> std::same_as<bool>;
    var.getTag();
  };

}  // namespace Ikarus::Concepts

#endif  // IKARUS_VARIABLE_INTERFACE_H
