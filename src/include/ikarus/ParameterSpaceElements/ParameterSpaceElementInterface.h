//
// Created by Alex on 27.05.2021.
//

#pragma once

namespace Ikarus::Concepts {
  template <typename ParameterSpaceElementType>
  concept ParameterSpaceElement = requires(ParameterSpaceElementType paraElement) {
    paraElement.dim();
    paraElement.numberOfEdges();
    paraElement.numberOfSurfaces();
    paraElement.numberOfVertices();
  };
}  // namespace Ikarus::Concepts