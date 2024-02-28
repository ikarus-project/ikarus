# Von Mises truss calculation

## Description

`iks007_vonMisesTruss.cpp` utilizes the tools and features mentioned in the previous examples to solve the
standard Von-Mises truss example found in literature (refer to Section 2[@misesTruss1923]).

## Code highlights

The struct named `Truss` is created such that it inherits from `FEBase`.
It must be decorated with `AutoDiffFE` as well to compute the stiffness matrix and load vectors during construction.
It is constructed as shown below:

```cpp
Truss(const BasisHandler &basisHandler, const typename LocalView::Element &element, double p_EA)
    : Base(basisHandler, element), EA{p_EA} {}
```

It takes a reference to the basis handler (`&basisHandler`), the element (`&element`), and the axial stiffness of the
truss structure (`p_EA`) as arguments during construction.

`ScalarType calculateScalarImpl(const FERequirementType &par, const Eigen::VectorX<ScalarType> &dx)` is
then defined, returning a scalar value, in this case the energy.
The energy is defined as `#!cpp 0.5 * EA * sqrt(LRefsquared) * Egl * Egl` with `Egl` being the Green-Lagrange strain defined as

```cpp
const Scalar Egl = 0.5 * (lsquared - LRefsquared) / LRefsquared;
```

The grid in this example is created explicitly. The grid is created here from the
[dune-foamgrid](https://www.dune-project.org/modules/dune-foamgrid/) module.
This allows the user to embed a one- or two-dimensional grid in a physical space of any dimension. Thus, here we embed
a one-dimensional truss system in a 2D plane. The height (`h`) and length (`L`) of the truss system are defined, which is followed by
the addition of the vertices and elements to create a grid as shown below:

```cpp
Dune::GridFactory<Dune::FoamGrid<1, 2, double>> gridFactory;
const double h = 1.0;
const double L = 2.0;
gridFactory.insertVertex({0, 0});
gridFactory.insertVertex({L, h});
gridFactory.insertVertex({2 * L, 0});
gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
auto grid     = gridFactory.createGrid();
auto gridView = grid->leafGridView();
```

The Lagrange basis is used to approximate the displacement field. The `Truss` elements are then created, followed by the fixing of the
degrees of freedom at the boundaries (`{0,0}` and `{2 * L,0}`). A vertical downward load is applied to the center node.
The non-linear operator is then constructed. The Newton-Raphson method is used as the non-linear solver, and an `nonLinearSolverObserver` is
created to write messages as desired by the non-linear solver. An additional `lvkObserver` is created using the `Ikarus::GenericControlObserver`
feature. This observer helps to fill up the matrix `lambdaAndDisp` with the load factor `lambda` and the two unconstrained degrees of
freedom whenever
the solution is changed (`#!cpp ControlMessages::SOLUTION_CHANGED`), which means that the Newton-Raphson method has converged to a solution.
This is implemented as depicted in the following:

```cpp
const int loadSteps = 10;
Eigen::Matrix3Xd lambdaAndDisp;
lambdaAndDisp.setZero(Eigen::NoChange, loadSteps + 1);
auto lvkObserver = std::make_shared<Ikarus::GenericObserver<Ikarus::ControlMessages>>(
    Ikarus::ControlMessages::SOLUTION_CHANGED, [&](int step) {
      lambdaAndDisp(0, step) = lambda;
      lambdaAndDisp(1, step) = d[2];
      lambdaAndDisp(2, step) = d[3];
    });
```

The load control method is used as the path-following strategy, and it is subscribed to both `vtkWriter` and `lvkObserver`.
The features from [Matplot++](https://github.com/alandefreitas/matplotplusplus) are then used to plot the load-displacement curve from
the matrix `lambdaAndDisp`.

## Takeaways

- `Dune::FoamGrid` can be used to embed one or two-dimensional grid entities into a multi-dimensional physical space.
- A simple truss element can be constructed using the automatic differentiation procedure.
- `Ikarus::GenericControlObserver` can be used to perform user-desired tasks at any desired point by observing a non-linear solver procedure.
