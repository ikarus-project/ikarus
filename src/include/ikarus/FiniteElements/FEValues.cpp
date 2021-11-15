//
// Created by lex on 14.11.21.
//

#include <ikarus/FiniteElements/FEValues.h>

namespace Ikarus::FiniteElements {
  const auto FEValues::get(const EntityType& eT, Ikarus::Variable::VariableTags varTag) const {
    auto tagFilter = [&varTag](auto&& var) { return isType(var, varTag); };
    switch (eT) {
      case EntityType::vertex:
        return vertexData | std::views::filter(tagFilter);
      case EntityType::edge:
        return edgeData | std::views::filter(tagFilter);
      case EntityType::surface:
        return surfaceData | std::views::filter(tagFilter);
      case EntityType::volume:
        return volumeData | std::views::filter(tagFilter);
    }
    __builtin_unreachable();
  }

  FEValues::EntityDataVectorType& FEValues::get(const EntityType& eT) {
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