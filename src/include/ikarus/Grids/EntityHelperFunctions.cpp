

#include <ikarus/Grids/EntityHelperFunctions.h>

namespace Ikarus {

  int determineEntityDimension(const EntityType& entType) {
    switch (entType) {
      case EntityType::vertex:
        return 0;
      case EntityType::edge:
        return 1;
      case EntityType::surface:
        return 2;
      case EntityType::volume:
        return 3;
      default:
        __builtin_unreachable();
    }
  }

  EntityType getEntityTypeFromCodim(int griddim, int codim) {
    assert(griddim - codim >= 0 && griddim - codim <= 3);
    switch (griddim - codim) {
      case 0:
        return Ikarus::EntityType::vertex;
      case 1:
        return Ikarus::EntityType::edge;
      case 2:
        return Ikarus::EntityType::surface;
      case 3:
        return Ikarus::EntityType::volume;
    }
    //    __builtin_unreachable();
    return EntityType::vertex;
  }
}  // namespace Ikarus