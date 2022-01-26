//
// Created by lex on 14.11.21.
//

#include <ikarus/FiniteElements/FEValues.h>

namespace Ikarus::FiniteElements {

  const FEValues::EntityDataVectorType& FEValues::get(const EntityType& eT) const {
    switch (eT) {
      case EntityType::vertex:
        return vertexData;
      case EntityType::edge:
        return edgeData;
      case EntityType::surface:
        return surfaceData;
      case EntityType::volume:
        return volumeData;
    }
    __builtin_unreachable();
  }

  void FEValues::add(const EntityType& eT, Ikarus::Variable::IVariable& var) {
    switch (eT) {
      case EntityType::vertex:
        vertexData.emplace_back(var);
        return;
      case EntityType::edge:
        edgeData.emplace_back(var);
        return;
      case EntityType::surface:
        surfaceData.emplace_back(var);
        return;
      case EntityType::volume:
        volumeData.emplace_back(var);
        return;
    }
    __builtin_unreachable();
  }

  void FEValues::add(const EntityType& eT, std::span<Ikarus::Variable::IVariable> varVec) {
    switch (eT) {
      case EntityType::vertex:
        vertexData.insert(vertexData.end(), varVec.begin(), varVec.end());
        return;
      case EntityType::edge:
        edgeData.insert(edgeData.end(), varVec.begin(), varVec.end());
        return;
      case EntityType::surface:
        surfaceData.insert(surfaceData.end(), varVec.begin(), varVec.end());
        return;
      case EntityType::volume:
        volumeData.insert(volumeData.end(), varVec.begin(), varVec.end());
        return;
    }
    __builtin_unreachable();
  }
}  // namespace Ikarus::FiniteElements