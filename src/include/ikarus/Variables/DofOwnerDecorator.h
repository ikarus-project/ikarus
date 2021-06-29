//
// Created by Alex on 16.06.2021.
//

#pragma once
#include <ikarus/Variables/GenericVariable.h>
#include <ikarus/Variables/VariableInterface.h>
#include <ikarus/utils/std/algorithms.h>
#include <ranges>

namespace Ikarus::Variable {
  class DofOwnerDecorator {
  public:
    void addDof(Ikarus::Concepts::Variable auto&& var) { Ikarus::stl::appendUnique(vars, Ikarus::Variable::GenericVariable(var)); }

    [[nodiscard]] auto& getDof(Ikarus::Concepts::Variable auto&& var) {
      // find using Tag
      auto it = std::ranges::find_if(
          vars, [](auto&& var1) { return (getTag(var1) == std::remove_cvref_t<decltype(var)>::tagvalue); });

      if (it == end(vars)) throw std::logic_error("Dof not available. Did you add it before with addDof?");

      return (*it);
    }

    std::vector<Ikarus::Variable::GenericVariable*> getDofs() {
      std::vector<Ikarus::Variable::GenericVariable*> pVars;

      for(auto&& var: vars)
        pVars.push_back(&var);

      return pVars;
    }

  private:
    std::vector<Ikarus::Variable::GenericVariable> vars;
  };

}  // namespace Ikarus::Variable