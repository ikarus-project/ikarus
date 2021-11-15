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
    const auto get(const EntityType& eT, Ikarus::Variable::VariableTags varTag) const;

    EntityDataVectorType& get(const EntityType& eT);

    void add(const EntityType& eT, Ikarus::Variable::IVariable& var);

    void add(const EntityType& eT, std::span<Ikarus::Variable::IVariable> varVec);

  private:
    EntityDataVectorType vertexData;
    EntityDataVectorType edgeData;
    EntityDataVectorType surfaceData;
    EntityDataVectorType volumeData;
  };
}  // namespace Ikarus::FiniteElements