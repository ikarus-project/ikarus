//
// Created by lex on 14.11.21.
//

#pragma once
#include <vector>

#include <ikarus/Grids/EntityHelperFunctions.h>
#include <ikarus/Variables/InterfaceVariable.h>

namespace Ikarus::FiniteElements {
  class FEValues {
    using EntityDataVectorType = std::vector<std::reference_wrapper<Ikarus::Variable::IVariable>>;

  public:
    auto get(const EntityType& eT, Ikarus::Variable::VariableTags varTag) const {
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

    const EntityDataVectorType& get(const EntityType& eT) const;

    void add(const EntityType& eT, Ikarus::Variable::IVariable& var);

    void add(const EntityType& eT, std::span<Ikarus::Variable::IVariable> varVec);

  private:
    EntityDataVectorType vertexData;
    EntityDataVectorType edgeData;
    EntityDataVectorType surfaceData;
    EntityDataVectorType volumeData;
  };
}  // namespace Ikarus::FiniteElements