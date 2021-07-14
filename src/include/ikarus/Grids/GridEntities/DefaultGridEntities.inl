
#include <ranges>
#include <span>

namespace Ikarus::Grid {

  template <int griddim, int wdim>
  [[nodiscard]] unsigned int DefaultGridEntity<griddim, 0, wdim>::subEntities(unsigned int codim) const {
    assert(codim <= 3 && codim > 0 && "Only subentities with 0< codim <= 3 supported");
    if constexpr (mydimension == 1)
      return getChildVertices().size();
    else if constexpr (mydimension == 2) {
      if (codim == 1)
        return getChildEntities<1>().size();
      else if (codim == 2)
        return getChildEntities<0>().size();
      else  // codim = 3
        throw std::logic_error("An entity with dimension 2 does not have subentities with codim 3!");
    } else {  //(mydimension == 3)
      if (codim == 1)
        return getChildEntities<2>().size();
      else if (codim == 2)
        return getChildEntities<1>().size();
      else  // codim = 3
        return getChildEntities<0>().size();
    }
  }

  template <int griddim, int cogriddim, int wdim>
  [[nodiscard]] unsigned int DefaultGridEntity<griddim, cogriddim, wdim>::subEntities(unsigned int codim) const {
    assert(codim > 0 && codim <= 2 && "Two dimensional entities only have subentities in 0<codimension<=2.");

    if constexpr (mydimension == 1) return getChildVertices().size();
    if constexpr (mydimension == 2) {
      if (codim == 1)
        return getChildEntities<1>().size();
      else if (codim == 2)
        return getChildEntities<0>().size();
    }
  }
  template <int griddim, int wdim>
  auto transformVertexPositionToDuneFieldVectorVector(std::span<DefaultGridEntity<griddim, 0, wdim>>&& vertices) {
    std::vector<Dune::FieldVector<double, wdim>> fieldVectorVector;
    for (const auto& pos :
         std::ranges::transform_view(vertices, &DefaultGridEntity<griddim, griddim, wdim>::getPosition))
      fieldVectorVector.push_back(toFieldVector(pos));
    return fieldVectorVector;
  }

  template <int griddim, int cogriddim, int wdim>
  auto DefaultGridEntity<griddim, cogriddim, wdim>::geometry() const {
    std::vector<Dune::FieldVector<double, dimensionworld>> fieldVectorVector
        = transformVertexPositionToDuneFieldVectorVector(getChildVertices());
    return Ikarus::Geometry::IGeometry(Geometry(type(), fieldVectorVector));
  }

  template <int griddim, int wdim>
  auto DefaultGridEntity<griddim, 0, wdim>::geometry() const {
    std::vector<Dune::FieldVector<double, dimensionworld>> fieldVectorVector
        = transformVertexPositionToDuneFieldVectorVector(getChildVertices());
    return Ikarus::Geometry::IGeometry(Geometry(type(), fieldVectorVector));
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(DefaultGridEntity<griddim, cogriddim, wdim>* const gridEntity) {
    return std::span(gridEntity->getChildVertices().begin(), gridEntity->getChildVertices().end());
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(DefaultGridEntity<griddim, cogriddim, wdim> const * const gridEntity) {
    return std::span(gridEntity->getChildVertices().begin(), gridEntity->getChildVertices().end());
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (cogriddim != griddim)
      return std::span(gridEntity.getChildVertices().begin(), gridEntity.getChildVertices().end());
    else
      throw std::logic_error("Vertices do not offer an iterator over vertices");
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(const DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (cogriddim != griddim)
      return std::span(gridEntity.getChildVertices().begin(), gridEntity.getChildVertices().end());
    else
      throw std::logic_error("Vertices do not offer an iterator over vertices");
  }

  template <int griddim, int cogriddim, int wdim>
  auto volumes(DefaultGridEntity<griddim, cogriddim, wdim>* gridEntity) {
    static_assert(cogriddim != 0, "Elements themself can not return span to iterate over themself");
    return std::span(gridEntity->getFatherElements().begin(), gridEntity->getFatherElements().end());
  }

  template <int griddim, int cogriddim, int wdim>
  auto volumes(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    static_assert(cogriddim != 0, "Elements themself can not return span to iterate over themself");
    return std::span(gridEntity.getFatherElements().begin(), gridEntity.getFatherElements().end());
  }

  template <int griddim, int cogriddim, int wdim>
  auto edges(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (griddim == cogriddim)  // gridEntity is a vertex!
      return std::span(gridEntity.template getFatherEntities<griddim - 1>().begin(),
                       gridEntity.template getFatherEntities<griddim - 1>().end());
    else
      return std::span(gridEntity.template getChildEntities<1>().begin(),
                       gridEntity.template getChildEntities<1>().end());
  }

  template <int griddim, int cogriddim, int wdim>
  auto edges(DefaultGridEntity<griddim, cogriddim, wdim>* gridEntity) {
    return edges(*gridEntity);
  }

  template <int griddim, int cogriddim, int wdim>
  auto surfaces(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (griddim == cogriddim)  // gridEntity is a vertex!
      return std::span(gridEntity.template getFatherEntities<griddim - 2>().begin(),
                       gridEntity.template getFatherEntities<griddim - 2>().end());
    else
      return std::span(gridEntity.template getChildEntities<2>().begin(),
                       gridEntity.template getChildEntities<2>().end());
  }

  template <int griddim, int cogriddim, int wdim, size_t dimE>
  requires requires { griddim >= dimE; }
  auto entities(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity, Dune::index_constant<dimE>&) {
    return std::span(gridEntity.template getChildEntities<dimE>().begin(),
                     gridEntity.template getChildEntities<dimE>().end());
  }

  /** \brief Return the fundamental geometric type of the entity, specialization for elements (codim==0) */
  template <int griddim, int wdim>
  Dune::GeometryType DefaultGridEntity<griddim, 0, wdim>::type() const {
    switch (mydimension) {
      case 0:
        DUNE_THROW(Dune::InvalidStateException, "The type of this element should not be only a vertex");
        ;
      case 1:
        return Dune::GeometryTypes::line;
      case 2:
        switch (subEntities(dimension)) {
          case 0:
            DUNE_THROW(Dune::InvalidStateException, "Type can not be deduced since verteces are not inserted.");
          case 2:
            DUNE_THROW(Dune::NotImplemented, "This should be a line and should be alread caught in case 1.");
          case 3:
            return Dune::GeometryTypes::triangle;
          case 4:
            return Dune::GeometryTypes::quadrilateral;
          default:
            DUNE_THROW(Dune::NotImplemented, "There is no fundamental geometry type for more than 4 vertices.");
        }
      case 3:
        switch (subEntities(dimension)) {
          case 4:
            return Dune::GeometryTypes::tetrahedron;
          case 5:
            return Dune::GeometryTypes::pyramid;
          case 6:
            return Dune::GeometryTypes::prism;
          case 8:
            return Dune::GeometryTypes::hexahedron;
          default:
            DUNE_THROW(Dune::NotImplemented, "There is no fundamental geometry type for these number of.");
        }
      default:
        DUNE_THROW(Dune::NotImplemented, "ERROR:  Unknown geometry type");
    }
  }

  /** \brief Return the fundamental geometric type of the entity, general  */
  template <int griddim, int cogriddim, int wdim>
  Dune::GeometryType DefaultGridEntity<griddim, cogriddim, wdim>::type() const {
    switch (mydimension) {
      case 0:
        DUNE_THROW(Dune::InvalidStateException, "The type of this subEntity should not be only a vertex");
        ;
      case 1:
        return Dune::GeometryTypes::line;
      case 2:
        switch (subEntities(dimension)) {
          case 0:
            DUNE_THROW(Dune::InvalidStateException, "Type can not be deduced since verteces are not inserted.");
          case 2:
            DUNE_THROW(Dune::NotImplemented, "This should be a line and should be alread caught in case 1.");
          case 3:
            return Dune::GeometryTypes::triangle;
          case 4:
            return Dune::GeometryTypes::quadrilateral;
          default:
            DUNE_THROW(Dune::NotImplemented, "There is no fundamental geometry type for more than 4 vertices.");
        }
      default:
        DUNE_THROW(Dune::NotImplemented, "ERROR:  Unknown geometry type");
    }
  }

}  // namespace Ikarus::Grid
