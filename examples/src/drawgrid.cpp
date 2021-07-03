
#include <dune/geometry/type.hh>

#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

int main() {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<3, 3>;
  SimpleGridFactory<Grid> gridFactory;
  using vertexType = Eigen::Vector3d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(vertexType{0.0, 0.0, -3.0});  // 0
  verticesVec.emplace_back(vertexType{2.0, 0.0, -3.0});  // 1
  verticesVec.emplace_back(vertexType{0.0, 2.0, -3.0});  // 2
  verticesVec.emplace_back(vertexType{2.0, 2.0, -3.0});  // 3
  verticesVec.emplace_back(vertexType{0.0, 0.0, 3.0});   // 4
  verticesVec.emplace_back(vertexType{2.0, 0.0, 3.0});   // 5
  verticesVec.emplace_back(vertexType{0.0, 2.0, 3.0});   // 6
  verticesVec.emplace_back(vertexType{2.0, 2.0, 3.0});   // 7
  verticesVec.emplace_back(vertexType{4.0, 0.0, 3.0});   // 8

  for (auto&& vert : verticesVec)
    gridFactory.insertVertex(vert);

  std::vector<size_t> elementIndices;
  elementIndices.resize(8);

  elementIndices = {0, 1, 2, 3, 4, 5, 6, 7};
  gridFactory.insertElement(Dune::GeometryTypes::hexahedron, elementIndices);
  elementIndices.resize(4);
  elementIndices = {1, 8, 3, 5};
  gridFactory.insertElement(Dune::GeometryTypes::tetrahedron, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();

  //  for(auto&& vertex : vertices(gridView))
  //    vertex.
  draw(gridView);
}