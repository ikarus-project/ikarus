
#include <dune/geometry/type.hh>

#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

int main() {
  using namespace Ikarus::Grid;
  SimpleGridFactory<3,3> gridFactory;
  using vertexType = Eigen::Vector3d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(0.0, 0.0, -3.0);  // 0
  verticesVec.emplace_back(2.0, 0.0, -3.0);  // 1
  verticesVec.emplace_back(0.0, 2.0, -3.0);  // 2
  verticesVec.emplace_back(2.0, 2.0, -3.0);  // 3
  verticesVec.emplace_back(0.0, 0.0, +3.0);  // 4
  verticesVec.emplace_back(2.0, 0.0, +3.0);  // 5
  verticesVec.emplace_back(0.0, 2.0, +3.0);  // 6
  verticesVec.emplace_back(2.0, 2.0, +3.0);  // 7
  verticesVec.emplace_back(4.0, 0.0, +3.0);  // 8

  for (auto&& vert : verticesVec)
    gridFactory.insertVertex(vert);

  std::vector<size_t> elementIndices;
  elementIndices.resize(8);

  elementIndices = {0, 1, 2, 3, 4, 5, 6, 7};
  gridFactory.insertElement(Ikarus::GeometryType::linearHexahedron, elementIndices);
  elementIndices.resize(4);
  elementIndices = {1, 8, 3, 5};
  gridFactory.insertElement(Ikarus::GeometryType::linearTetrahedron, elementIndices);

  auto grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();

  draw(gridView);
}