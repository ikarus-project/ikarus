// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/vtk/datacollectors/unstructureddatacollector.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/lagrangepoints.hh>
#include <dune/geometry/virtualrefinement.hh>
namespace Dune::Vtk {
  /// Implementation of \ref Discontinuous DataCollector for Iga cells with trimming information
  template <class GridView>
  requires(GridView::dimension == 2 and GridView::dimensionworld==3) class Shell2DDataCollector
    : public UnstructuredDataCollectorInterface<GridView, Shell2DDataCollector<GridView>, Partitions::All> {
    using Self  = Shell2DDataCollector;
    using Super = UnstructuredDataCollectorInterface<GridView, Self, Partitions::All>;


    constexpr static int dimw = GridView::dimensionworld;

    typedef typename GridView::Grid::ctype ctype;
  public:
    using Super::dim;
    typedef VirtualRefinement<2, ctype> Refinement;
    typedef typename Refinement::IndexVector IndexVector;
    typedef typename Refinement::ElementIterator SubElementIterator;
    typedef typename Refinement::VertexIterator SubVertexIterator;




    using Super::partition;
    Dune::RefinementIntervals subSampleInPlane;



  public:
    explicit Shell2DDataCollector(GridView const& gridView, Dune::RefinementIntervals subSampleInPlane_ = RefinementIntervals(0)) :
   Super(gridView),subSampleInPlane{subSampleInPlane_} {
    }


    auto determineRefinementType(const Dune::GeometryType& type) const
    {
      if(type==Dune::GeometryTypes::quadrilateral)
        return Dune::GeometryTypes::quadrilateral;
      else if(type==Dune::GeometryTypes::triangle)
        return Dune::GeometryTypes::triangle;
      else
        DUNE_THROW(Dune::NotImplemented,"Shell2DDataCollector works only for element types that are quads or triangles");
    }


    /// Construct the point sets
    void updateImpl() {
      pointSets_.clear();
      numPoints_ = 0;
      numCells_  = 0;
      std::int64_t vertex_idx = 0;
      auto const& indexSet    = gridView_.indexSet();
      for (auto const& eit : elements(gridView_)) {
        const auto refElementType= determineRefinementType(eit.type());
        Refinement& inPlaneRefinement = buildRefinement<2, ctype>(refElementType, refElementType);

        for (SubVertexIterator inPlaneVertexIt = inPlaneRefinement.vBegin(subSampleInPlane), send = inPlaneRefinement.vEnd(subSampleInPlane);
             inPlaneVertexIt != send; ++inPlaneVertexIt) {
            vertexIndex_.emplace(std::array<int, 2>({indexSet.index(eit), inPlaneVertexIt.index()}), vertex_idx++);
            ++numPoints_;
          }
          numCells_+=inPlaneRefinement.nElements(subSampleInPlane);
        }
    }

    /// Return number of Lagrange nodes
    [[nodiscard]] std::uint64_t numPointsImpl() const { return numPoints_; }

    /// Return a vector of point coordinates.
    /**
     * The vector of point coordinates is composed of vertex coordinates of the untrimmed elements and
     * the vertices of the triangulated trimmed elements
     **/
    template <class T>
    [[nodiscard]] std::vector<T> pointsImpl() const {
      std::vector<T> data(this->numPoints() * 3);
      auto const& indexSet = gridView_.indexSet();
      for (auto eit : elements(gridView_)) {
        auto geometry          = eit.geometry();
        const int elementId = indexSet.index(eit);

        const auto refElementType= determineRefinementType(eit.type());
        Refinement& refinement = buildRefinement<2, ctype>(refElementType, refElementType);
        for (SubVertexIterator sit = refinement.vBegin(subSampleInPlane), send = refinement.vEnd(subSampleInPlane);
             sit != send; ++sit) {
         auto v          = geometry.global(sit.coords());
          std::int64_t idx = 3 * vertexIndex_.at({elementId,sit.index()});
          for (std::size_t j = 0; j < v.size(); ++j)
            data[idx + j] = T(v[j]);
          for (std::size_t j = v.size(); j < 3u; ++j)
            data[idx + j] = T(0);
        }
      }

      return data;
    }

    /// Return number of grid cells
    [[nodiscard]] std::uint64_t numCellsImpl() const { return numCells_; }


    //this should be a frre function
    void reorder(const std::span<int64_t>& vA, const std::span<size_t>& vI) const
    {
      size_t i, j, k;
      for(i = 0; i < vA.size(); i++){
        while(i != (j = vI[i])){
          k = vI[j];
          std::swap(vA[j], vA[k]);
          std::swap(vI[i], vI[j]);
        }
      }
    }
    /// \brief Return cell types, offsets, and connectivity. \see Cells
    /**
     * The cell connectivity is composed of cell vertices
     **/
    [[nodiscard]] Cells cellsImpl() const {
      Cells cells;
      cells.connectivity.reserve(this->numPoints());
      cells.offsets.reserve(this->numCells());
      cells.types.reserve(this->numCells());

      auto const& indexSet = gridView_.indexSet();

      std::int64_t old_o = 0;
      for (auto const& eit : elements(gridView_, partition)) {
        const int elementId = indexSet.index(eit);
        const auto refElementType= determineRefinementType(eit.type());

        Refinement& refinement = buildRefinement<2, ctype>(refElementType, refElementType);
        for (SubElementIterator eRit = refinement.eBegin(subSampleInPlane), eRend = refinement.eEnd(subSampleInPlane);
             eRit != eRend; ++eRit) {
          Vtk::CellType cellType(refElementType, Vtk::CellType::LINEAR);

          const size_t numberOfVirtualVerticesInElement = eRit.vertexIndices().size();
            for (auto idxLocal : eRit.vertexIndices()) {
              std::int64_t idx
                  = vertexIndex_.at({elementId, idxLocal});
              cells.connectivity.push_back(idx);
            }
            std::vector<size_t> permIndices;
            for (size_t i = 0; i < numberOfVirtualVerticesInElement; ++i) {
              permIndices.push_back(cellType.permutation(i));
            }

              auto cellConBegin= cells.connectivity.end()-numberOfVirtualVerticesInElement;
              reorder({cellConBegin,cells.connectivity.end()},permIndices);

          cells.types.push_back(cellType.type());
          cells.offsets.push_back(old_o += numberOfVirtualVerticesInElement);
        }
      }
      return cells;
    }

    /// Evaluate the `fct` at element vertices and edge centers in the same order as the point coords.
    template <class T, class GlobalFunction>
    [[nodiscard]] std::vector<T> pointDataImpl(GlobalFunction const& fct) const {
      int nComps = fct.numComponents();
      std::vector<T> data(this->numPoints() * nComps);
            auto localFct        = localFunction(fct);
      auto const& indexSet = gridView_.indexSet();
      for (auto eit : elements(gridView_)) {
        const int elementId = indexSet.index(eit);
        localFct.bind(eit);
        const auto refElementType= determineRefinementType(eit.type());
        Refinement& refinement = buildRefinement<2, ctype>(refElementType, refElementType);
        for (SubVertexIterator sit = refinement.vBegin(subSampleInPlane), send = refinement.vEnd(subSampleInPlane);
             sit != send; ++sit) {

          std::int64_t idx = nComps * vertexIndex_.at({elementId,sit.index()});

          for (std::size_t comp = 0; comp < nComps; ++comp)
              data[idx + comp] = T(localFct.evaluate(comp, sit.coords()));
        }
      }

      return data;
    }
    // Evaluate `fct` in center of cell.
    template <class T, class VtkFunction>
    [[nodiscard]] std::vector<T> cellDataImpl(VtkFunction const& fct) const {
      int nComps = fct.numComponents();
      std::vector<T> data;
      data.reserve(this->numCells_ * nComps);
//
      auto localFct        = localFunction(fct);
      auto const& indexSet = gridView_.indexSet();
      for (auto eit : elements(gridView_)) {
        localFct.bind(eit);

        auto geometry          = eit.geometry();
        const int elementId = indexSet.index(eit);

        const auto refElementType= determineRefinementType(eit.type());
        Refinement& refinement = buildRefinement<2, ctype>(refElementType, refElementType);
        for (SubElementIterator eVit = refinement.eBegin(subSampleInPlane), eend = refinement.eEnd(subSampleInPlane);
             eVit != eend; ++eVit) {

          for (std::size_t comp = 0; comp < nComps; ++comp)
              data.push_back(T(localFct.evaluate(comp, eVit.coords())));
        }
      }
      return data;
    }

  protected:
    using Super::gridView_;

    std::uint64_t numPoints_ = 0;
    std::uint64_t numCells_  = 0;

    using PointSet = LagrangePointSet<typename GridView::ctype, GridView::dimension>;
    std::map<GeometryType, PointSet> pointSets_;
//    std::vector<std::int64_t> indexMap_;
    std::map<std::array<int , 2>, std::int64_t> vertexIndex_;
  };

}  // namespace Dune::Vtk
