//
// Created by Alex on 26.05.2021.
//

#pragma once
#include <ranges>
#include <string>
#include <vector>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>

#include <ikarus/Geometries/GeometryInterface.h>
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>
#include <ikarus/utils/utils/traits.h>

namespace Ikarus::Grid {

  template <int dim, int dimworld>
  class SimpleGrid;

  template <int griddim, int cogriddim, int wdim>
  class DefaultGridEntity;

  namespace Impl {
    /** \brief  The ChildEntityPointerTupleGenerator generates a tuple of grid subentities, which can be of different
     * size depending on the dimension of the entity, e.g. a entity of dimension 3 has surfaces, edges and vertices as
     * children whereas a entity with dimension 1 only has vertices as children
     */
    template <int griddim, int mydim, int wdim, int... codim>
    static std::tuple<std::vector<DefaultGridEntity<griddim, griddim - codim, wdim>*>...>
        ChildEntityPointerTupleGenerator(std::integer_sequence<int, codim...>);

    /** \brief  The FatherEntityPointerTupleGenerator generates a tuple of grid entities, which can be of different
     * size depending on the codimension of the entity, e.g. a entity with codimension 3 has volumes, surfaces and edges
     * as fathers whereas a entity with codimension 1 only has edges as fathers
     */
    template <int griddim, int mydim, int wdim, int... codim>
    static std::tuple<std::vector<DefaultGridEntity<griddim, codim, wdim>*>...> FatherEntityPointerTupleGenerator(
        std::integer_sequence<int, codim...>);
  }  // namespace Impl

  /**
   * \brief DefaultGridEntity
   *
   * \tparam  griddim  The dimension of the grid
   * \tparam  cogriddim=0  The codimension of the entity
   * \tparam  wdim  The dimension of the world space where the grid is embedded
   *
   * \note Partial template specialization for entities with codim == 0
   * These entities have no father
   **/
  template <int griddim, int wdim>
  class DefaultGridEntity<griddim, 0, wdim> {
  public:
    static constexpr int dimension      = griddim;
    static constexpr int codimension    = 0;
    static constexpr int mydimension    = dimension;
    static constexpr int dimensionworld = wdim;

    DefaultGridEntity(int levelInput, size_t idInput) : levelIndex{levelInput}, id{idInput} {}

    auto& getChildVertices() { return std::get<0>(entitiesChildren); }
    const auto& getChildVertices() const { return std::get<0>(entitiesChildren); }

    template <int dimEnt>
    auto& getChildEntities() {
      static_assert(dimEnt >= 0 && dimEnt < griddim, "You asked for a non-existing ChildEntity!");
      return std::get<dimEnt>(entitiesChildren);
    }

    template <int dimEnt>
    const auto& getChildEntities() const {
      static_assert(dimEnt >= 0 && dimEnt < griddim, "You asked for a non-existing ChildEntity!");
      return std::get<dimEnt>(entitiesChildren);
    }

    /** \brief Type of the geometry of this entity */
    using Geometry = Dune::MultiLinearGeometry<double, codimension, dimensionworld>;

    /** \brief Type of the containter for the grid childrens of this entity, it stores only pointers to the entities. */
    using EntitiesChildernType
        = decltype(Impl::ChildEntityPointerTupleGenerator<dimension, mydimension, dimensionworld>(
            std::make_integer_sequence<int, mydimension>()));

    /** \brief Returns the number of subEntities of this entity, e.g. a line has two vertices as
     * subtypes */
    [[nodiscard]] unsigned int subEntities(unsigned int codim) const;

    /** \brief Return the fundamental geometric type of the entity */
    [[nodiscard]] Ikarus::GeometryType type() const;

    /** \brief Returns the geometric realization of the entity */
    auto geometry() const;

    /** \brief Get refínement level where this entity belongs to*/
    [[nodiscard]] int level() const { return levelIndex; }

  private:
    friend class SimpleGrid<griddim, wdim>;
    /** \brief The refinement level to which this entity belongs */
    int levelIndex{};
    /** \brief A persistent id of this entity*/
    size_t id{};

    /** \brief Return copy of the id of this entity */
    [[nodiscard]] size_t getID() const { return id; }

    /** \brief Childrens of the entity on the current grid , i.e. surfaces of a cube*/
    EntitiesChildernType entitiesChildren;

    /** \brief Childrens of the entity on a finer grid , i.e. subcubes of a cube, it stores only pointers to the
     * entities.*/
    std::vector<DefaultGridEntity<dimension, 0, dimensionworld>*> levelEntitiesChildren;
  };

  /**
   * \brief DefaultGridEntity
   *
   * \tparam  griddim  The dimension of the grid
   * \tparam  cogriddim=griddim  The codimension of the entity. Here
   * \tparam  wdim  The dimension of the world space where the grid is embedded
   *
   * \note Partial template specialization for entities with codim == griddim
   * These entities have no grid children, since they are vertices!
   **/
  template <int griddim, int wdim>
  class DefaultGridEntity<griddim, griddim, wdim> {
  public:
    static constexpr int dimension      = griddim;
    static constexpr int codimension    = griddim;
    static constexpr int mydimension    = 0;
    static constexpr int dimensionworld = wdim;

    DefaultGridEntity(int levelInput, const Eigen::Vector<double, wdim>& vecInput, size_t idInput)
        : levelIndex{levelInput}, id{idInput}, position{vecInput} {}

    auto& getFatherElements() { return std::get<0>(entitiesFathers); }

    template <int coDimsubEnt>
    auto& getFatherEntities() {
      static_assert(coDimsubEnt >= 0 && coDimsubEnt <= griddim - 1,
                    "You asked for a non-existing FatherEntity!");  // TODO Fatherentities
      return std::get<coDimsubEnt>(entitiesFathers);
    }
    DefaultGridEntity() = default;

    /** \brief Type of the geometry of this entity */
    using Geometry = Dune::MultiLinearGeometry<double, codimension, dimensionworld>;

    /** \brief Type of the containter for the grid fathers of this entity */
    using EntitiesFatherType = decltype(Impl::FatherEntityPointerTupleGenerator<dimension, mydimension, dimensionworld>(
        std::make_integer_sequence<int, codimension>()));

    /** \brief Get refínement level where this entity belongs to*/
    [[nodiscard]] int level() const { return levelIndex; }

    /** \brief Return position of this vertex */
    const Eigen::Vector<double, wdim>& getPosition() { return position; }

    /** \brief Return the fundamental geometric type of the entity */
    [[nodiscard]] Ikarus::GeometryType type() const { return Ikarus::GeometryType::vertex; }

    /** \brief Returns the number of subEntities of this entity, e.g. a line has two verteces as
     * subtypes */
    [[nodiscard]] unsigned int subEntities(unsigned int) const { return 0; }

    /** \brief Returns the geometric realization of the entity */
    auto geometry() const { return Geometry(duneType(type()), position); }

  private:
    friend class SimpleGrid<griddim, wdim>;
    /** \brief The index of this element on the level it belongs to */
    int levelIndex{};
    /** \brief A persistent id of this entity*/
    size_t id{};

    /** \brief Return copy of the id of this entity */
    [[nodiscard]] size_t getID() const { return id; }

    /** \brief Childrens of the entity on the current grid , i.e. surfaces of a cube*/
    EntitiesFatherType entitiesFathers;

    /** \brief Childrens of the entity on a finer grid , i.e. subcubes of a cube
     * Raw pointers to the quantities belonging to the grid are stored here */
    std::vector<DefaultGridEntity<dimension, codimension, dimensionworld>*> levelEntitiesChildren;

    /** \brief The position of the vertex */
    Eigen::Vector<double, wdim> position{};
  };

  /**
   * \brief DefaultGridEntity
   *
   * \note Partial template specialization for entities with codim != griddim and codim != 0
   * These entities live inbetween the vertices and elements, i.e. edges,surface
   **/
  template <int griddim, int cogriddim, int wdim>
  class DefaultGridEntity {
  public:
    static constexpr int dimension      = griddim;
    static constexpr int codimension    = cogriddim;
    static constexpr int mydimension    = griddim - cogriddim;
    static constexpr int dimensionworld = wdim;

    DefaultGridEntity(int levelInput, size_t idInput) : levelIndex{levelInput}, id{idInput} {}

    /** \brief Get refinement level where this entity belongs to*/
    [[nodiscard]] int level() const { return levelIndex; }

    /** \brief Type of the geometry of this entity */
    using Geometry = Dune::MultiLinearGeometry<double, codimension, dimensionworld>;

    template <int dimEnt>
    auto& getChildEntities() {
      static_assert(dimEnt >= 0 && dimEnt < griddim, "You asked for a non-existing ChildEntity!");
      return std::get<dimEnt>(entitiesChildren);
    }

    template <int dimEnt>
    const auto& getChildEntities() const {
      static_assert(dimEnt >= 0 && dimEnt < griddim, "You asked for a non-existing ChildEntity!");
      return std::get<dimEnt>(entitiesChildren);
    }

    const auto& getChildVertices() const { return std::get<0>(entitiesChildren); }
    auto& getChildVertices() { return std::get<0>(entitiesChildren); }

    /** \brief Returns the number of subEntities of this entity, e.g. a line has two verteces as
     * subtypes */
    [[nodiscard]] unsigned int subEntities(unsigned int codim) const;

    /** \brief Return the fundamental geometric type of the entity */
    [[nodiscard]] Ikarus::GeometryType type() const;

    /** \brief Returns the geometric realization of the entity */
    auto geometry() const;

    /** \brief Return copy of the id of this entity */
  private:
    friend class SimpleGrid<griddim, wdim>;
    /** \brief The refinement level to which this entity belongs */
    int levelIndex{};

    /** \brief A persistent id of this entity*/
    size_t id{};

    /** \brief Return copy of the id of this entity */
    [[nodiscard]] size_t getID() const { return id; }

    /** \brief Childrens of the entity on the current grid , i.e. surfaces of a cube*/
    decltype(Impl::ChildEntityPointerTupleGenerator<dimension, mydimension, dimensionworld>(
        std::make_integer_sequence<int, mydimension>())) entitiesChildren;

    /** \brief Childrens of the entity on a finer grid , i.e. subcubes of a cube*/
    std::vector<DefaultGridEntity<dimension, codimension, dimensionworld>*> levelEntitiesChildren;
  };
}  // namespace Ikarus::Grid

#include "DefaultGridEntities.inl"