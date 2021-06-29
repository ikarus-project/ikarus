
#include <ranges>
#include <span>

namespace Ikarus::Grid {

  template <int griddim, int cogriddim, int wdim>
  auto vertices(DefaultGridEntity<griddim, cogriddim, wdim>* gridEntity) {
    return std::span(gridEntity->getChildVertices().begin(), gridEntity->getChildVertices().end());
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(DefaultGridEntity<griddim, cogriddim, wdim> const * const  gridEntity) {
    return std::span(gridEntity->getChildVertices().begin(), gridEntity->getChildVertices().end());
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (cogriddim!=griddim)
    return std::span(gridEntity.getChildVertices().begin(), gridEntity.getChildVertices().end());
    else
      throw std::logic_error("Vertices do not offer an iterator over vertices");
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(const DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (cogriddim!=griddim)
      return std::span(gridEntity.getChildVertices().begin(), gridEntity.getChildVertices().end());
    else
      throw std::logic_error("Vertices do not offer an iterator over vertices");
  }

  // TODO: Check if possible With Children and Father inversed

  template <int griddim, int cogriddim, int wdim>
  auto elements(DefaultGridEntity<griddim, cogriddim, wdim>* gridEntity) {
    static_assert(cogriddim != 0, "Elements themself can not return span to iterate over themself");
    return std::span(gridEntity->getFatherElements().begin(), gridEntity->getFatherElements().end());
  }

  template <int griddim, int cogriddim, int wdim> auto edges(DefaultGridEntity<griddim, cogriddim, wdim>* gridEntity) {
    if constexpr (griddim == cogriddim)  // gridEntity is a vertex!
      return std::span(gridEntity->template getFatherEntities<griddim - 1>().begin(),
                       gridEntity->template getFatherEntities<griddim - 1>().end());
    else
      return std::span(gridEntity->template getChildEntities<1>().begin(),
                       gridEntity->template getChildEntities<1>().end());
  }

  template <int griddim, int cogriddim, int wdim> auto edges(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (griddim == cogriddim)  // gridEntity is a vertex!
      return std::span(gridEntity.template getFatherEntities<griddim - 1>().begin(),
                       gridEntity.template getFatherEntities<griddim - 1>().end());
    else
      return std::span(gridEntity.template getChildEntities<1>().begin(),
                       gridEntity.template getChildEntities<1>().end());
  }

  template <int griddim, int cogriddim, int wdim>
  auto surfaces(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    return gridEntity.template getChildEntities<1>();
  }

  /** \brief Return the fundamental geometric type of the entity, specialization for elements (codim==0) */
  template <int griddim, int wdim> Dune::GeometryType DefaultGridEntity<griddim, 0, wdim>::type() const {
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
