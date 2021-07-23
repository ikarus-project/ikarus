//
// Created by Alex on 16.06.2021.
//

#pragma once
#include <ranges>

#include <ikarus/Variables/InterfaceVariable.h>
#include <ikarus/Variables/VariablePolicies.h>
#include <ikarus/utils/utils/algorithms.h>

namespace Ikarus::Variable {
  class DofOwnerDecorator {
  public:
    void addDof(Ikarus::Concepts::Variable auto&& var) {
      Ikarus::utils::appendUnique(vars, Ikarus::Variable::IVariable(var));
    }

    [[nodiscard]] auto& getDof(Ikarus::Concepts::Variable auto&& var) {
      // find using Tag
      auto it = std::ranges::find_if(
          vars, [](auto&& var1) { return (getTag(var1) == std::remove_cvref_t<decltype(var)>::tagvalue); });

      if (it == end(vars)) throw std::logic_error("Dof not available. Did you add it before with addDof?");

      return (*it);
    }

    std::vector<Ikarus::Variable::IVariable*> getDofs() {
      std::vector<Ikarus::Variable::IVariable*> pVars;

      for (auto&& var : vars)
        pVars.push_back(&var);

      return pVars;
    }

  private:
    std::vector<Ikarus::Variable::IVariable> vars;
  };

}  // namespace Ikarus::Variable