//
// Created by Alex on 28.09.2021.
//

#pragma once
#include <map>

#include <ikarus/Grids/GridViews/SimpleGridView.h>

template <typename GridType>
class SimpleGridIndexSet {
public:
  static constexpr int dimension      = GridType::dimension;
  static constexpr int worldDimension = GridType::dimensionworld;

  SimpleGridIndexSet(GridType& grid, int level) : grid_{&grid}, level_{level} {
    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<dimension + 1>()), [&](auto i) {
      size_t indexPerEntityOfCodim = 0;
      for (auto& entity : grid.gridEntitiesContainer.template getSubEntities<i>(level)) {
        indexContainer_[i][grid.getEntityID(entity)] = indexPerEntityOfCodim++;
      }
    });
  }

  template <typename Entity>
  size_t index(const Entity& entity) const {
    return indexContainer_[std::remove_pointer_t<Entity>::codimension].at(grid_->getEntityID(entity));
  }

  template <typename Entity>
  size_t subIndex(const Entity& entity, int subEntityNumber, unsigned int codim) const {
    constexpr size_t entityDimension = std::remove_pointer_t<Entity>::mydimension;
    //    constexpr size_t worldDimension  = std::remove_pointer_t<Entity>::dimensionworld;
    assert(codim <= dimension);
    if constexpr (entityDimension > 1)
      if (dimension - codim == 1) return index(*entity.template getChildEntities<1>()[subEntityNumber]);  // edges

    if constexpr (entityDimension > 2)
      if (dimension - codim == 2) return index(*entity.template getChildEntities<2>()[subEntityNumber]);  // surfaces

    if (dimension - codim == 0)
      return index(*entity.template getChildEntities<0>()[subEntityNumber]);  // these are vertices
    else                                                                      // codim==0
    {
      assert(subEntityNumber == 0 && "Your choice of entity and codim means that you request the id the entity itself. This is only valid for i==0 since each entity has only one unique index");
      return index(entity);
    }
  }

  int size(int codim) const { return indexContainer_[codim].size(); }

private:
  // in the first index the rootEntities are stored
  std::array<std::map<size_t, size_t>, dimension + 1> indexContainer_;
  GridType* grid_;
  int level_;
};
