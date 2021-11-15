//
// Created by lex on 05.10.21.
//

#pragma once
#include <optional>
#include <ranges>
#include <vector>

#include <ikarus/Variables/InterfaceVariable.h>

namespace Ikarus {
  enum class EntityType : int { vertex, edge, surface, volume };
  template <int subEntity>
  struct EntitiesWithCoDim {};

  int determineEntityDimension(const EntityType& entType);

  EntityType getEntityTypeFromCodim(int griddim, int codim);
}  // namespace Ikarus