//
// Created by lex on 05.10.21.
//

#pragma once
#include <optional>
#include <ranges>
#include <vector>

#include <ikarus/Variables/InterfaceVariable.h>

namespace Ikarus {
  enum class EntityType : int { vertex, edge, surface, volume, generic };

  template <int subEntity>
  struct EntitiesWithCoDim {};
  struct VerticesStruct {};
  struct EdgesStruct {};
  struct SurfacesStruct {};
  struct VolumesStruct {};
  struct RootEntitiesStruct {};

  inline constexpr auto Vertices     = VerticesStruct{};
  inline constexpr auto Edges        = EdgesStruct{};
  inline constexpr auto Surfaces     = SurfacesStruct{};
  inline constexpr auto Volumes      = VolumesStruct{};
  inline constexpr auto RootEntities = RootEntitiesStruct{};

  template <typename EntityTypeG, int gridDim>
  static consteval int determineEntityDimension() {
    using EntityType = std::decay_t<EntityTypeG>;
    //  std::cout<<Dune::className(EntityType{})<<std::endl;
    if constexpr (std::is_same_v<EntityType, VerticesStruct> || std::is_same_v<EntityType, EntitiesWithCoDim<gridDim>>)
      return 0;
    else if constexpr (std::is_same_v<EntityType,
                                      EdgesStruct> || std::is_same_v<EntityType, EntitiesWithCoDim<gridDim - 1>>)
      return 1;
    else if constexpr (std::is_same_v<EntityType,
                                      SurfacesStruct> || std::is_same_v<EntityType, EntitiesWithCoDim<gridDim - 2>>)
      return 2;
    else if constexpr (std::is_same_v<EntityType,
                                      VolumesStruct> || std::is_same_v<EntityType, EntitiesWithCoDim<gridDim - 3>>)
      return 3;
    else  // if constexpr (std::is_same_v<EntityType, RootEntitiesStruct> || std::is_same_v<EntityType,
          // EntitiesWithCoDim<0>>)
      return gridDim;
  }

  class FEValues {
    std::vector<std::reference_wrapper<Ikarus::Variable::IVariable>> vertexData;
    std::vector<std::reference_wrapper<Ikarus::Variable::IVariable>> edgeData;
    std::vector<std::reference_wrapper<Ikarus::Variable::IVariable>> surfaceData;
    std::vector<std::reference_wrapper<Ikarus::Variable::IVariable>> volumeData;
    std::vector<std::reference_wrapper<Ikarus::Variable::IVariable>> genericData;

  public:
    const auto get(const EntityType& eT, Ikarus::Variable::VariablesTags varTag) const {
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
        case EntityType::generic:
          return genericData | std::views::filter(tagFilter);
      }
    }

    auto& get(const EntityType& eT) {
      switch (eT) {
        case EntityType::vertex:
          return vertexData;
        case EntityType::edge:
          return edgeData;
        case EntityType::surface:
          return surfaceData;
        case EntityType::volume:
          return volumeData;
        default:  //     EntityType::generic:
          return genericData;
      }
    }

    void set(const EntityType& eT, Ikarus::Variable::IVariable& var) {
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
        case EntityType::generic:
          genericData.emplace_back(var);
          return;
      }
    }

    void set(const EntityType& eT, std::span<Ikarus::Variable::IVariable> varVec) {
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
        case EntityType::generic:
          genericData.insert(genericData.end(), varVec.begin(), varVec.end());
          return;
      }
    }
  };

}  // namespace Ikarus