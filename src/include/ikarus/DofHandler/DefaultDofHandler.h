//
// Created by Alex on 26.05.2021.
//

#pragma once

#include <unordered_set>

#include <ikarus/FiniteElements/PhysicalElementPolicies.h>
#include <ikarus/Variables/VariableDefinitions.h>

#include <ikarus/utils/std/hashs.h>
namespace Ikarus::Variable {
  class GenericVariable;
}

namespace Ikarus::DofHandler {

  template <typename FEContainer>
  requires Concepts::HasgetDofVector<typename std::remove_cvref_t<FEContainer>::value_type> || Concepts::
      HasFreegetDofVector<typename std::remove_cvref_t<FEContainer>::value_type>
  class DefaultDofHandler {
  public:
    using DofSet      = std::unordered_set<Ikarus::Variable::GenericVariable*>;
    using DofVector   = std::vector<Ikarus::Variable::GenericVariable*>;

    DefaultDofHandler(FEContainer& feContainer) : feContainer_{feContainer} {}

    void createDofList()  {
      DofSet dofSet;
      for (auto&& fe : feContainer_)
      {
        auto eleDofs = getDofVector(fe);
        dofSet.insert(begin(eleDofs), end(eleDofs));
      }

      for(auto&& dof: dofSet) {
//     auto pair = std::make_pair(Ikarus::Grid::IGridEntity(dof.first),dof.second)
        dofVec.push_back(dof);
      }
    }

  private:
    DofVector dofVec;
    FEContainer&  feContainer_;
  };
}  // namespace Ikarus::DofHandler