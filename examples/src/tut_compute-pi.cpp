//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>
#define ALUGRIDDEBUG false
#include <dune/alugrid/grid.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/boundarysegment.hh>
//#include <dune/common/function.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/Grids/GridHelper/griddrawer.h>

/**
 * \brief Mapping from R^d to the graph of a given function
 */
template <int dim, int dimworld>
class GraphMapping {
  // Element corners in the global parameter domain
  std::array<Dune::FieldVector<double, dim>, dim + 1> corners_;

  std::function<Dune::FieldVector<double, 2>(Dune::FieldVector<double, 2>)> graph_;

public:
  GraphMapping(std::array<Dune::FieldVector<double, dim>, dim + 1> corners,
               std::function<Dune::FieldVector<double, 2>(Dune::FieldVector<double, 2>)> graph)
      : corners_(corners), graph_(graph) {}

  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  Dune::FieldVector<double, dimworld> operator()(const Dune::FieldVector<double, dim>& x) const {
    // Linear interpolation between the corners
    auto globalX = corners_[0];
    for (size_t i = 0; i < x.size(); i++)
      for (int j = 0; j < dim; j++)
        globalX[j] += x[i] * (corners_[i + 1][j] - corners_[0][j]);
    return graph_(globalX);
  }
};

struct BS : Dune::BoundarySegment<2, 2, double> {

  BS(const Dune::FieldVector<double, 2>& a ,const Dune::FieldVector<double, 2>& b) : corners{{a,b}} {}
  virtual Dune::FieldVector<double, 2> operator()(const Dune::FieldVector<double, 1>& local) const override {
    double angle =  + local[0] * (2 * std::numbers::pi - 0);

    Dune::FieldVector<double, 2> result = {0, 0};
    result = corners[0] + local[0]*(corners[1]-corners[0]);

    return result/result.two_norm();
  }

  std::array<Dune::FieldVector<double, 2>, 2> corners;
};

int main(int argc, char** argv) {
  // use MPI helper to initialize MPI
  Dune::MPIHelper::instance(argc, argv);

  using namespace Dune;
  constexpr int gridDim = 2;
  Dune::GridFactory<Dune::ALUGrid<gridDim, 2, Dune::cube, Dune::nonconforming>> gridFactory;
  const double h = 1.0;
  const double L = 1.0;
  std::array<FieldVector<double, 2>, 4> corners0
      = {{{-sqrt(2) / 2, -sqrt(2) / 2}, {sqrt(2) / 2, -sqrt(2) / 2}, {sqrt(2) / 2, sqrt(2) / 2}, {-sqrt(2) / 2, sqrt(2) / 2}}};
  gridFactory.insertVertex(corners0[0]);
  gridFactory.insertVertex(corners0[1]);
  gridFactory.insertVertex(corners0[2]);
  gridFactory.insertVertex(corners0[3]);

  auto parametrization = [](const FieldVector<double, 2>& x) -> FieldVector<double, 2> { return x / x.two_norm(); };

  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1, 3,2});
  gridFactory.insertBoundarySegment({0, 1}, std::make_shared<BS>(corners0[0],corners0[1]));
  gridFactory.insertBoundarySegment({1, 2}, std::make_shared<BS>(corners0[1],corners0[2]));
  gridFactory.insertBoundarySegment({0, 3}, std::make_shared<BS>(corners0[0],corners0[3]));
  gridFactory.insertBoundarySegment({3, 2}, std::make_shared<BS>(corners0[3],corners0[2]));


  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();
  draw(gridView);

  double area1 = 0.0;
  for (auto element : elements(gridView)) {
    area1 += element.geometry().volume();
  }

  auto f       = [](auto&& global) { return sqrt(global[0] * global[0] + global[1] * global[1]); };
  double area2 = 0.0;
  for (auto& element : elements(gridView)) {
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(element.type(), 10, Dune::QuadratureType::GaussLegendre);
    for (auto& gp : rule)
      area2 += element.geometry().integrationElement(gp.position()) * gp.weight();
    //      area2 += f(element.geometry().global(gp.position())) * element.geometry().integrationElement(gp.position())
    //               * gp.weight();  // integrationElement --> JacobiDeterminant
  }

  std::cout << "Area2 " << area2 << " " << std::numbers::pi << std::endl;
  bool adaptive = true;
  for (int i = 0; i < 10; ++i) {
    area1 = 0.0;
    grid->globalRefine(1);

    auto gridViewRefined = grid->leafGridView();
    std::cout << "This gridview contains: ";
    std::cout << gridViewRefined.size(0) << " elements" << std::endl;

//    draw(gridViewRefined);
    for (auto& element : elements(gridViewRefined)) {
      area1 += element.geometry().volume();
    }
    std::cout << "area1 " << area1 << " " << std::numbers::pi << std::endl;
    draw(gridViewRefined);
  }

  double circumference = 0.0;
  for (auto& element : elements(grid->leafGridView()))
    if (element.hasBoundaryIntersections())
      for (auto& intersection : intersections(grid->leafGridView(), element))
        if (intersection.boundary()) circumference += intersection.geometry().volume();

  std::cout << circumference << std::endl;

}