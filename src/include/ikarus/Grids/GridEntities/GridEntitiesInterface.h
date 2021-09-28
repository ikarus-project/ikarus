//
// Created by Alex on 26.05.2021.
//

#pragma once
#include <ranges>
#include <span>

#include <dune/geometry/type.hh>

#include <ikarus/Geometries/GeometryInterface.h>

#include <ikarus/FiniteElements/FiniteElementInterface.h>
#include <ikarus/FiniteElements/PhysicalElementPolicies.h>
#include <ikarus/Geometries/GeometryInterface.h>

namespace Ikarus::Grid {

  /** \brief A type-erased grid Entity
   * This can be an element, a vertex or something else*/
  class IGridEntity {
  public:
    template <Ikarus::Concepts::GridEntity GE>
    explicit IGridEntity(const GE& fe) : gEimpl{std::make_unique<GEImpl<GE> >(fe)} {}

    ~IGridEntity() = default;
    IGridEntity(const IGridEntity& other) : gEimpl{other.gEimpl->clone()} {}
    IGridEntity& operator=(const IGridEntity& other) {
      IGridEntity tmp(other);  // Temporary-swap idiom
      std::swap(gEimpl, tmp.gEimpl);
      return *this;
    }

    IGridEntity(IGridEntity&&) noexcept = default;
    IGridEntity& operator=(IGridEntity&&) noexcept = default;

  private:
    struct GEBase {
      virtual ~GEBase()                                                     = default;
      [[nodiscard]] virtual int do_level() const                            = 0;
      [[nodiscard]] virtual size_t do_getID() const                         = 0;
      [[nodiscard]] virtual Ikarus::Geometry::IGeometry do_geometry() const = 0;
      [[nodiscard]] virtual Dune::GeometryType do_type() const              = 0;
      [[nodiscard]] virtual std::span<IGridEntity> do_vertices() const      = 0;
      [[nodiscard]] virtual std::span<IGridEntity> do_edges() const         = 0;
      [[nodiscard]] virtual std::span<IGridEntity> do_vertices()            = 0;
      [[nodiscard]] virtual std::span<IGridEntity> do_edges()               = 0;
      [[nodiscard]] virtual std::unique_ptr<GEBase> clone() const           = 0;  // Prototype Design Pattern
    };

    template <typename GE>
    struct GEImpl : public GEBase {
      template <typename GridEntity>
      auto transformSpan(GridEntity* fe) {
        auto createIGE = [](GridEntity* ge) { return IGridEntity(ge); };
        return std::ranges::transform_view(vertices(fe), createIGE);
      }

      template <typename GridEntity>
      auto transformSpan(GridEntity& fe) {
        auto createIGE = [](GridEntity& ge) { return IGridEntity(ge); };
        return std::ranges::transform_view(vertices(fe), createIGE);
      }

      template <typename GridEntity>
      auto transformSpan(const GridEntity& fe) const {
        auto createIGE = [](const GridEntity& ge) { return IGridEntity(ge); };
        return std::ranges::transform_view(vertices(fe), createIGE);
      }

      explicit GEImpl(GE gearg) : gE{gearg} {};
      [[nodiscard]] int do_level() const final { return gE.level(); }
      [[nodiscard]] size_t do_getID() const final { return gE.getID(); }
      [[nodiscard]] Ikarus::Geometry::IGeometry do_geometry() const final { return gE.geometry(); }
      [[nodiscard]] Dune::GeometryType do_type() const final { return gE.type(); };
      [[nodiscard]] std::span<IGridEntity> do_vertices() const final { return transformSpan(gE); };
      [[nodiscard]] std::span<IGridEntity> do_edges() const final { return edges(gE); };
      [[nodiscard]] std::span<IGridEntity> do_vertices() final { return transformSpan(gE); };
      [[nodiscard]] std::span<IGridEntity> do_edges() final { return edges(gE); };
      [[nodiscard]] std::unique_ptr<GEBase> clone() const final { return std::make_unique<GEImpl>(*this); }
      GE gE;
    };

    std::unique_ptr<GEBase> gEimpl;  // Pimpl idiom / Bridge Design Patterns

    friend int level(const IGridEntity& fe);
    friend auto geometry(const IGridEntity& fe);
    friend auto type(const IGridEntity& fe);
    friend auto vertices(const IGridEntity& fe);
    friend std::span<const IGridEntity> edges(const IGridEntity& fe);
    friend auto vertices(IGridEntity& fe);
    friend std::span<IGridEntity> edges(IGridEntity& fe);
    friend std::span<const IGridEntity> vertices(const IGridEntity* fe);
    friend std::span<const IGridEntity> edges(const IGridEntity* fe);
    friend std::span<IGridEntity> vertices(IGridEntity* fe);
    friend std::span<IGridEntity> edges(IGridEntity* fe);
    friend auto getID(const IGridEntity& fe);
  };

  int level(const IGridEntity& fe) { return fe.gEimpl->do_level(); }
  auto geometry(const IGridEntity& fe) { return fe.gEimpl->do_geometry(); }
  auto type(const IGridEntity& fe) { return fe.gEimpl->do_type(); }
  auto vertices(const IGridEntity& fe) { return fe.gEimpl->do_vertices(); }
  std::span<const IGridEntity> edges(const IGridEntity& fe) { return fe.gEimpl->do_edges(); }
  auto vertices(IGridEntity& fe) {
    return fe.gEimpl->do_vertices();
    ;
  }
  std::span<IGridEntity> edges(IGridEntity& fe) { return fe.gEimpl->do_edges(); }
  std::span<const IGridEntity> vertices(const IGridEntity* fe) { return fe->gEimpl->do_vertices(); }
  std::span<const IGridEntity> edges(const IGridEntity* fe) { return fe->gEimpl->do_edges(); }
  std::span<IGridEntity> vertices(IGridEntity* fe) { return fe->gEimpl->do_vertices(); }
  std::span<IGridEntity> edges(IGridEntity* fe) { return fe->gEimpl->do_edges(); }
  auto getID(const IGridEntity& fe) { return fe.gEimpl->do_getID(); }

}  // namespace Ikarus::Grid
