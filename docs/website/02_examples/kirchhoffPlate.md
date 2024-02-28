# Plate subjected to a surface load

## Description

Kirchhoff-type plate element is implemented in `iks004_kirchhoffPlate.cpp` using the automatic differentiation
technique as commented before. The basis used for discretization is a NURBS basis from the `dune-iga` module.
The problem is solved, and convergence plots are created by comparing the solutions to the available analytical solutions for the
simply supported case.

## Code highlights

Similar to the `struct` named `Solid` in `iks003_incompressible_LinearElasticity.cpp`, here a `struct` named `KirchhoffPlate` is created.
It inherits from `FEBase` and it must be decorated with `AutoDiffFE` as well to compute the required matrices and vectors.
It is constructed as shown below:

```cpp
KirchhoffPlate(const BasisHandler &basisHandler, const typename LocalView::Element &element, double p_Emodul,
               double p_nu, double p_thickness)
    : Base(basisHandler, element),
      Emodul{p_Emodul},
      nu{p_nu},
      thickness{p_thickness} {
  geometry_.emplace(this->localView().element().geometry());
}
```

It takes in the `p_thickness` parameter in addition to the ones in `Solid`. Here, the energy is calculated as:

```cpp
energy += (0.5 * kappa.dot(D * kappa) - w * lambda) * geometry_->integrationElement(gp.position()) * gp.weight();
```

with `kappa` being the vector of curvature containing $\kappa_{xx}, \kappa_{yy}$ and $\kappa_{xy}$. If the boundaries
are clamped, the penalty method could be used, for example, to fix the derivatives of the displacements. Here, however, the example for
the simply supported case is discussed further.

A two-dimensional NURBS grid is created from the [dune-iga](https://github.com/rath3t/dune-iga) module.

```cpp
constexpr int griddim                                    = 2; // (1)!
constexpr int dimworld                                   = 2; // (2)!
const std::array<std::vector<double>, griddim> knotSpans = { { {0, 0, 1, 1}, {0, 0, 1, 1} } }; // (3)!
using ControlPoint = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointType;
const double Lx = 10; // (4)!
const double Ly = 10; // (5)!

const std::vector<std::vector<ControlPoint>> controlPoints
    = { { {.p = {0, 0}, .w = 1}, {.p = {0, Ly}, .w = 1} },
        { {.p = {Lx, 0}, .w = 1}, {.p = {Lx, Ly}, .w = 1} } }; // (6)!

std::array<int, griddim> dimsize = {2, 2}; // (7)!

auto controlNet = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointNetType(dimsize, controlPoints); // (8)!
using Grid      = Dune::IGA::NURBSGrid<griddim, dimworld>;

Dune::IGA::NURBSPatchData<griddim, dimworld> patchData;

patchData.knotSpans     = knotSpans; // (9)!
patchData.degree        = {1, 1}; // (10)!
patchData.controlPoints = controlNet; // (11)!
patchData = Dune::IGA::degreeElevate(patchData, 0, 1); // (12)!
patchData = Dune::IGA::degreeElevate(patchData, 1, 1); // (13)!
Grid grid(patchData); // (14)!
```

1. The dimension of the grid. Here, two-dimensional.
2. The dimension of the physical space in which the grid lies (the embedding space of the grid). Here, two-dimensional.
3. The knot vector for the NURBS grid.
4. Length of the plate
5. Width of the plate
6. Control points and the weights for the square plate (the control points are ordered such that all the control points for a particular
   $x$-position are listed first, followed by the subsequent $x$-positions in ascending order).
7. Number of control points in either direction.
8. Creation of the control net.
9. Binding the `knotSpans` to a particular NURBS `patchData`.
10. The polynomial degree for the basis in either direction. It can also be calculated from the `knotSpans`.
11. Binding the `controlNet` to a particular NURBS `patchData`.
12. Elevating the polynomial degree for the `patchData` in $x$-direction (`0`) by 1.
13. Elevating the polynomial degree for the `patchData` in $y$-direction (`1`) by 1.
14. Creating the grid object from the patch data

In order to obtain the convergence plots, the system is solved five times, with the refinement level
increasing by 1 each time using the command `#!cpp grid.globalRefine(1);`.
The NURBS basis can be obtained from the freestanding functions `nurbs()`, as shown below:

```cpp
auto basis = Ikarus::makeBasis(gridView, nurbs());
```

This is followed by specifying the Dirichlet boundary conditions, creating the finite elements and the assembler,
solving the system of equations, and post-processing using Paraview as mentioned in the previous examples.
The analytical solution for the simply supported case is adapted from
[Wikipedia](https://en.wikipedia.org/wiki/Bending_of_plates#Simply-supported_plate_with_uniformly-distributed_load) and is also mentioned below:

```cpp
auto wAna = [&](auto x) {
  double w                = 0.0;
  const int seriesFactors = 40;
  const double pi         = std::numbers::pi;
  auto oddFactors
      = std::ranges::iota_view(1, seriesFactors) | std::views::filter([](auto i) { return i % 2 != 0; });
  for (auto m : oddFactors)
    for (auto n : oddFactors)
      w += sin(m * pi * x[0] / Lx) * sin(n * pi * x[1] / Ly)
           / (m * n * Dune::power(m * m / (Lx * Lx) + n * n / (Ly * Ly), 2));

  return 16 * totalLoad / (Dune::power(pi, 6) * D) * w;
};
```

The `#!cpp Dune::Functions::makeDiscreteGlobalBasisFunction` is used to create a function from the nodal finite element
solution of the displacements and the NURBS basis whereas the `#!cpp Dune::Functions::makeAnalyticGridViewFunction` is
used to create a function by using the function to evaluate the analytical solutions and the `gridView` to get the position `x`.
Local functions are then created that are used later to calculate the $L^2$-error.

```cpp
auto wGlobalFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 1>>(basis.flat(), w);
auto wGlobalAnalyticFunction = Dune::Functions::makeAnalyticGridViewFunction(wAna, gridView);
auto localw                  = localFunction(wGlobalFunction);
auto localwAna               = localFunction(wGlobalAnalyticFunction);
```

The $L^2$-error is calculated by using
$$
L^2\textrm{-error} = \frac{\sqrt{\sum_{ele} \int_{\Omega_{ele}} \left( w_{analytical}-w_{FE} \right)^2}}{L^2\textrm{-exact}}
$$
with $L^2\textrm{-exact} = \sqrt{\sum_{ele} \int_{\Omega_{ele}} \left( w_{analytical}\right)^2}$
as shown below:

```cpp
double l2_error = 0.0;
double l2_normEx = 0.0;
for (auto &ele : elements(gridView)) {
  localView.bind(ele);
  localw.bind(ele);
  localwAna.bind(ele);
  const auto geo   = localView.element().geometry();
  const auto &rule = Dune::QuadratureRules<double, 2>::rule(
      ele.type(), 2U * localView.tree().finiteElement().localBasis().order());
  for (auto gp : rule) {
    const auto intElement = ele.geometry().integrationElement(gp.position()) * gp.weight();
    const auto w_ex       = localwAna(gp.position());
    const auto w_fe       = localw(gp.position());
    l2_error += Dune::power(w_ex - w_fe, 2) * intElement;
    l2_normEx += w_ex * intElement;
  }
}
l2_error = std::sqrt(l2_error) / std::sqrt(l2_normEx);
```

The number of degrees of freedom for each refinement level and its corresponding $L^2$-error is pushed to a vector that
can be later used to create plots using the features from Matlab.

```cpp
std::vector<size_t> dofsVec;
std::vector<double> l2Evcector;
dofsVec.push_back(basis.flat().size());
l2Evcector.push_back(l2_error);
```

## Takeaways

- NURBS grids can be created using the `dune-iga` module.
- The basis for the corresponding NURBS grid can be obtained using the `nurbs()` function.
- The Kirchhoff plate element can be easily implemented by evaluating the energy and using automatic differentiation methods.
- $L^2$-error can be evaluated to perform convergence studies.
