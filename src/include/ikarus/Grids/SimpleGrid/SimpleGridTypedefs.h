//
// Created by ac120950 on 20.07.2021.
//

#pragma once

namespace SimpleGridTypedefs{
  template<int dimensionworld>
  struct VertexIndexPair {
    Eigen::Vector<double, dimensionworld> vertex;
    size_t index;
  };
}
