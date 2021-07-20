//
// Created by ac120950 on 20.07.2021.
//

#ifndef IKARUS_SIMPLEGRIDTYPEDEFS_H
#define IKARUS_SIMPLEGRIDTYPEDEFS_H

namespace SimpleGridTypedefs{
  template<int dimensionworld>
  struct VertexIndexPair {
    Eigen::Vector<double, dimensionworld> vertex;
    size_t index;
  };
}



#endif  // IKARUS_SIMPLEGRIDTYPEDEFS_H
