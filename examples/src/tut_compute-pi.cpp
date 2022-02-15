//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/geometry/quadraturerules.hh>
//#include <dune/common/function.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/Grids/GridHelper/griddrawer.h>

/**
 * \brief Mapping from R^d to the graph of a given function
 */
template<int dim, int dimworld>
class GraphMapping
{
  // Element corners in the global parameter domain
  std::array<Dune::FieldVector<double,dim>, dim+1> corners_;

  std::function<Dune::FieldVector<double,2>(Dune::FieldVector<double,2>)> graph_;

public:
  GraphMapping(std::array<Dune::FieldVector<double,dim>, dim+1> corners, std::function<Dune::FieldVector<double,2>(Dune::FieldVector<double,2 > ) > graph)
      : corners_(corners), graph_(graph)
  {}

  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  Dune::FieldVector<double,dimworld> operator() (const Dune::FieldVector<double,dim>& x) const
  {
    // Linear interpolation between the corners
    auto globalX = corners_[0];
    for (size_t i=0; i<x.size(); i++)
      for (int j=0; j<dim; j++)
        globalX[j] += x[i]*(corners_[i+1][j]-corners_[0][j]);
    return graph_(globalX);
  }

};

int main() {
  using namespace Dune;
  constexpr int gridDim = 2;
  Dune::GridFactory<Dune::FoamGrid<2, 2, double>> gridFactory;
  const double h = 1.0;
  const double L = 1.0;
  std::array<FieldVector<double,2>, 3> corners0 = {{{-sqrt(2)/2, -sqrt(2)/2}, {sqrt(2)/2, -sqrt(2)/2}, {0, 1}}};
  gridFactory.insertVertex(corners0[0]);
  gridFactory.insertVertex(corners0[1]);
  gridFactory.insertVertex(corners0[2]);

  auto parametrization = [](const FieldVector<double,2>& x) -> FieldVector<double,2>
  {return x/x.two_norm();};

  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0,1,2}, GraphMapping<gridDim, gridDim>(corners0, parametrization));

  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();
  draw(gridView);

  double area1 = 0.0;
  for (auto& element : elements(gridView)) {
    area1 += element.geometry().volume();
  }

  auto f       = [](auto&& global) { return sqrt(global[0] * global[0] + global[1] * global[1]); };
  double area2 = 0.0;
  for (auto& element : elements(gridView)) {
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(element.type(), 1, Dune::QuadratureType::GaussLegendre);
    for (auto& gp : rule)
      //      area2 += element.geometry().integrationElement(gp.position()) * gp.weight();
      area2 += f(element.geometry().global(gp.position())) * element.geometry().integrationElement(gp.position())
               * gp.weight();  // integrationElement --> JacobiDeterminant
  }

  std::cout << area1 << " " << area2 << std::endl;

  for (int i = 0; i < 3; ++i) {
    area1 = 0.0;
    grid->globalRefine(1);

    auto gridViewRefined = grid->leafGridView();
    std::cout << "This gridview contains: ";
    std::cout << gridViewRefined.size(0) << " elements" << std::endl;

    draw(gridViewRefined);
    for (auto& element : elements(gridViewRefined)) {
      area1 += element.geometry().volume();
    }
    std::cout << area1 << " " << std::numbers::pi << std::endl;
  }

  double circumference = 0.0;
  for (auto& element : elements(gridView))
    if (element.hasBoundaryIntersections())
      for (auto& intersection : intersections(gridView, element))
        if (intersection.boundary()) circumference += intersection.geometry().volume();

  std::cout << circumference << std::endl;
}