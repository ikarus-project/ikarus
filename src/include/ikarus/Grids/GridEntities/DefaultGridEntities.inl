
#include <ranges>
#include <span>

#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/utils/utils/algorithms.h>

namespace Ikarus::Grid {

  template <int griddim, int wdim>
  [[nodiscard]] unsigned int DefaultGridEntity<griddim, 0, wdim>::subEntities(unsigned int codim) const {
    if (codim > 3 || codim == 0) throw std::logic_error("Only subentities with 0< codim <= 3 supported");
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
    assert(codim > 1 && codim <= 3 && "Two dimensional entities only have subentities in 0<codimension<=2.");

    if constexpr (mydimension == 1) return getChildVertices().size();
    if constexpr (mydimension == 2) {
      if (wdim - codim == 1)
        return getChildEntities<1>().size();
      else if (wdim - codim == 0)
        return getChildEntities<0>().size();
    }
    return 0;
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
    return Ikarus::Geometry::IGeometry(Geometry(duneType(type()), fieldVectorVector));
  }

  template <int griddim, int wdim>
  auto DefaultGridEntity<griddim, 0, wdim>::geometry() const {
    std::vector<Dune::FieldVector<double, dimensionworld>> fieldVectorVector
        = transformVertexPositionToDuneFieldVectorVector(getChildVertices());
    return Ikarus::Geometry::IGeometry(Geometry(duneType(type()), fieldVectorVector));
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (cogriddim != griddim)
      return Ikarus::utils::transformPointerRangeToReferenceRange(gridEntity.getChildVertices());
    else
      throw std::logic_error("Vertices do not offer an iterator over vertices");
  }

  template <int griddim, int cogriddim, int wdim>
  auto vertices(const DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (cogriddim != griddim)
      return Ikarus::utils::transformPointerRangeToReferenceRange(gridEntity.getChildVertices());
    else
      throw std::logic_error("Vertices do not offer an iterator over vertices");
  }

  template <int griddim, int cogriddim, int wdim>
  auto volumes(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    static_assert(cogriddim != 0, "Elements themself can not return span to iterate over themself");
    return Ikarus::utils::transformPointerRangeToReferenceRange(gridEntity.getFatherElements());
  }

  template <int griddim, int cogriddim, int wdim>
  auto edges(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (griddim == cogriddim)  // gridEntity is a vertex!
      return Ikarus::utils::transformPointerRangeToReferenceRange(gridEntity.template getFatherEntities<griddim - 1>());
    else
      return Ikarus::utils::transformPointerRangeToReferenceRange(gridEntity.template getChildEntities<1>());
  }

  template <int griddim, int cogriddim, int wdim>
  auto surfaces(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity) {
    if constexpr (griddim == cogriddim)  // gridEntity is a vertex!
      return Ikarus::utils::transformPointerRangeToReferenceRange(gridEntity.template getFatherEntities<griddim - 2>());
    else
      return Ikarus::utils::transformPointerRangeToReferenceRange(gridEntity.template getChildEntities<2>());
  }

  template <int griddim, int cogriddim, int wdim, size_t dimE>
  requires requires { griddim >= dimE; }
  auto entities(DefaultGridEntity<griddim, cogriddim, wdim>& gridEntity, Dune::index_constant<dimE>&) {
    return Ikarus::utils::transformPointerRangeToReferenceRange(gridEntity.template getChildEntities<dimE>());
  }

  /** \brief Return the fundamental geometric type of the entity, specialization for elements (codim==0) */
  template <int griddim, int wdim>
  Ikarus::GeometryType DefaultGridEntity<griddim, 0, wdim>::type() const {
    switch (mydimension) {
      case 0:
        DUNE_THROW(Dune::InvalidStateException, "The type of this element should not be only a vertex");
        ;
      case 1:
        return Ikarus::GeometryType::linearLine;
      case 2:
        switch (subEntities(dimension)) {
          case 0:
            DUNE_THROW(Dune::InvalidStateException, "Type can not be deduced since verteces are not inserted.");
          case 2:
            DUNE_THROW(Dune::NotImplemented, "This should be a line and should be alread caught in case 1.");
          case 3:
            return Ikarus::GeometryType::linearTriangle;
          case 4:
            return Ikarus::GeometryType::linearQuadrilateral;
          default:
            DUNE_THROW(Dune::NotImplemented, "There is no fundamental geometry type for more than 4 vertices.");
        }
      case 3:
        switch (subEntities(dimension)) {
          case 4:
            return Ikarus::GeometryType::linearTetrahedron;
          case 5:
            return Ikarus::GeometryType::pyramid;
          case 6:
            return Ikarus::GeometryType::prism;
          case 8:
            return Ikarus::GeometryType::linearHexahedron;
          default:
            DUNE_THROW(Dune::NotImplemented, "There is no fundamental geometry type for these number of.");
        }
      default:
        DUNE_THROW(Dune::NotImplemented, "ERROR:  Unknown geometry type");
    }
  }

  /** \brief Return the fundamental geometric type of the entity, general  */
  template <int griddim, int cogriddim, int wdim>
  Ikarus::GeometryType DefaultGridEntity<griddim, cogriddim, wdim>::type() const {
    switch (mydimension) {
      case 0:
        DUNE_THROW(Dune::InvalidStateException, "The type of this subEntity should not be only a vertex");
        ;
      case 1:
        return Ikarus::GeometryType::linearLine;
      case 2:
        switch (subEntities(dimension)) {
          case 0:
            DUNE_THROW(Dune::InvalidStateException, "Type can not be deduced since verteces are not inserted.");
          case 2:
            DUNE_THROW(Dune::NotImplemented, "This should be a line and should be alread caught in case 1.");
          case 3:
            return Ikarus::GeometryType::linearTriangle;
          case 4:
            return Ikarus::GeometryType::linearQuadrilateral;
          default:
            DUNE_THROW(Dune::NotImplemented, "There is no fundamental geometry type for more than 4 vertices.");
        }
      default:
        DUNE_THROW(Dune::NotImplemented, "ERROR:  Unknown geometry type");
    }
  }
}  // namespace Ikarus::Grid
