//
// Created by ac120950 on 17.07.2020.
//

#pragma once

#include <Eigen/Core>
// namespace Ikarus::Concepts {
//  template <class GeometryEntityType>
//  concept Geometry = requires(
//      GeometryEntityType&& geoEntity,
//      Eigen::Matrix<typename GeometryEntityType::ctype, Eigen::Dynamic, GeometryEntityType::mydimension> dN,
//      Eigen::Matrix<typename GeometryEntityType::ctype, GeometryEntityType::worlddim, Eigen::Dynamic> x) {
////    { geoEntity.determinantJacobian(dN, x) } -> std::same_as<typename GeometryEntityType::ctype>;
////    { geoEntity.jacobianTransposed(dN, x) } -> std::same_as<typename GeometryEntityType::JacobianTransposed>;
////    {
////      geoEntity.jacobianInverseTransposed(dN, x)
////      } -> std::same_as<typename GeometryEntityType::JacobianInverseTransposed>;
//  };
//}  // namespace Ikarus::Concepts

namespace Ikarus::Geometry {

  /** \brief A type-erased grid Entity
   * This can be an element, a vertex or something else*/
  class IGeometry {
  public:
    template <typename GE>
    explicit IGeometry(const GE& fe) : gEimpl{std::make_unique<GEImpl<GE> >(fe)} {}

    ~IGeometry() = default;
    IGeometry(const IGeometry& other) : gEimpl{other.gEimpl->clone()} {}
    IGeometry& operator=(const IGeometry& other) {
      IGeometry tmp(other);  // Temporary-swap idiom
      std::swap(gEimpl, tmp.gEimpl);
      return *this;
    }

    IGeometry(IGeometry&&) noexcept = default;
    IGeometry& operator=(IGeometry&&) noexcept = default;

  private:
    struct GEBase {
      virtual ~GEBase() = default;
      //      virtual void do_initialize() = 0;
      //      [[nodiscard]] virtual int do_level() const = 0;
      //      [[nodiscard]] virtual Ikarus::Geometry::IGeometry do_geometry() const = 0;
      //      [[nodiscard]] virtual Dune::GeometryType do_type() const = 0;
      //      [[nodiscard]] virtual std::span<IGeometry> do_vertices() const = 0;
      //      [[nodiscard]] virtual std::span<IGeometry> do_edges() const = 0;
      [[nodiscard]] virtual std::unique_ptr<GEBase> clone() const = 0;  // Prototype Design Pattern
    };

    template <typename GE>
    struct GEImpl : public GEBase {
      explicit GEImpl(GE gearg) : gE{gearg} {};
      //      [[nodiscard]] int do_level() const final { return gE.level(); }
      //      [[nodiscard]] Ikarus::Geometry::IGeometry do_geometry() const final { return gE.geometry(); }
      //      [[nodiscard]] Dune::GeometryType do_type() const final { return gE.type(); };
      //      [[nodiscard]] std::span<IGeometry> do_vertices() const final { return vertices(gE); };
      //      [[nodiscard]] std::span<IGeometry> do_edges() const final { return edges(gE); };
      [[nodiscard]] std::unique_ptr<GEBase> clone() const final { return std::make_unique<GEImpl>(*this); }
      GE gE;
    };

    std::unique_ptr<GEBase> gEimpl;  // Pimpl idiom / Bridge Design Patterns

    //    friend int level(const IGeometry& fe);
    //    friend auto geometry(const IGeometry& fe);
    //    friend auto type(const IGeometry& fe);
    //    friend auto vertices(const IGeometry& fe);
    //    friend auto edges(const IGeometry& fe);
  };

  //  int level(const IGeometry& fe) { return fe.gEimpl->do_level(); }
  //  auto geometry(const IGeometry& fe) { return fe.gEimpl->do_geometry(); }
  //  auto type(const IGeometry& fe) { return fe.gEimpl->do_type(); }
  //  auto vertices(const IGeometry& fe) { return fe.gEimpl->do_vertices(); }
  //  auto edges(const IGeometry& fe) { return fe.gEimpl->do_vertices(); }

}  // namespace Ikarus::Geometry
