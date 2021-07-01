//
// Created by Alex on 19.05.2021.
//

#pragma once
#include <ostream>

namespace Ikarus::Concepts {
  template <typename ManifoldType>
  concept Manifold = requires(ManifoldType var, typename ManifoldType::CorrectionType correction, std::ostream& s,
                              typename ManifoldType::CoordinateType value) {
    typename ManifoldType::ctype;
    ManifoldType::valueSize;
    ManifoldType::correctionSize;
    typename ManifoldType::CoordinateType;
    typename ManifoldType::CorrectionType;
    { var.getValue() } -> std::convertible_to<typename ManifoldType::CoordinateType>;
    { var.setValue(value) } -> std::same_as<void>;
    { var.update(correction) } -> std::same_as<void>;
    { s << var } -> std::same_as<std::ostream&>;
  };
}  // namespace Ikarus::Concepts