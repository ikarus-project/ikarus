# Deformation of an incompressible rubber block

## Description

`iks003_incompressible_LinearElasticity.cpp` uses finite element technology with displacement and pressure as
independent degrees of freedom to simulate the deformation of an incompressible rubber block. The potential energy
for such a system is defined in the `Solid struct` by the function
`calculateScalarImpl(const FERequirementType &par, const Eigen::VectorX<ScalarType> &dx)`.
This function uses the principles of automatic differentiation to provide the stiffness matrix and
other necessary quantities to perform a static structural analysis.
It inherits from `FEBase` which provides information about the `localView` of the element.

## Code highlights

The `struct` named `Solid` is created. It is constructed as shown below:

```cpp
Solid(const BasisHandler &basisHandler, const typename LocalView::Element &element, double emod, double nu)
      : Base(basisHandler, element), emod_{emod}, nu_{nu} {
    mu_       = emod_ / (2 * (1 + nu_));
    lambdaMat = convertLameConstants({.emodul = emod_, .nu = nu_}).toLamesFirstParameter();
}
```

It takes a reference to the basis handler (`&basisHandler`),
the element (`&element`), and the material parameters, namely Young's modulus
(`emod`) and Poisson's ratio (`nu`), as arguments during construction.
The function `convertLameConstants()` is a helper function
to switch between the Lame parameters.

`ScalarType calculateScalarImpl(const FERequirementType &par, const Eigen::VectorX<ScalarType> &dx)` is
then defined, returning a scalar value, in this case the energy.
The energy is then calculated as follows:

```cpp
energy += (0.5 * (2 * mu_ * symgradu.squaredNorm() - 1 / lambdaMat * Dune::power(pressure, 2)) + pressure * divU
           - x.dot(fext))
          * geo.integrationElement(gp.position()) * gp.weight();  // plane strain for 2D
```

Here:

- `symgradu` is the symmetric part of the gradient of displacements
- `lambdaMat` is the first Lame parameter
- `pressure` and `x` are the nodal pressure and current position, respectively
- `divU` is the divergence of the displacement vector
- `fext` is the external force vector
- `gp.position()` and `gp.weight()` are the positions and weights from the quadrature rule
- `geo.integrationElement()` returns the determinant of Jacobian required from the iso-parametric concept

A Yasp 2D grid[@sander2020dune] is created of the size $L$ x $h$ with 20 elements in both directions, as shown below:

```cpp
using Grid        = Dune::YaspGrid<gridDim>;
const double L    = 1;
const double h    = 1;
const size_t elex = 20;
const size_t eley = 20;

Dune::FieldVector<double, 2> bbox       = {L, h};
std::array<int, 2> elementsPerDirection = {elex, eley};
auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);
auto gridView                           = grid->leafGridView();
```

A linear Lagrangian basis is opted for the displacements and a constant basis for the pressure degrees of freedom using the
`composite` basis feature from Dune, as shown below:

```cpp
auto basis = Ikarus::makeBasis(
gridView, composite(power<2>(lagrange<1>()), lagrange<0>()));
```

Here, `#!cpp power<2>` is used to approximate the displacement field in both $x$ and $y$ directions.
A vector of `Solid` finite elements that are decorated by `AutoDiffFE` are then constructed as shown below:

```cpp
std::vector<AutoDiffFE<Solid<decltype(basis)>>> fes;
for (auto &ele : elements(gridView))
  fes.emplace_back(basis, ele, Emod, nu);
```

The displacement degrees of freedom at position $y=0$ are fixed using the following snippet:

```cpp
auto basisP = std::make_shared<const decltype(basis)>(basis);
Ikarus::DirichletValues dirichletValues(basisP->flat());
dirichletValues.fixDOFs([](auto &basis_, auto &dirichletFlags) {
  Dune::Functions::forEachBoundaryDOF(subspaceBasis(basis_, _0),
                   [&](auto &&localIndex, auto &&localView, auto &&intersection) {
                       if (std::abs(intersection.geometry().center()[1]) < 1e-8)
                         dirichletFlags[localView.index(localIndex)] = true;
  });
});
```

Here, all the element edges lying on the boundary of the physical domain are looped over and checked to see if the first index
of the center of the edge (`intersection.geometry().center()[1]`) is close to zero. If this is the case, the corresponding $x$-displacement
degrees of freedom (obtained via `subspaceBasis(basis_, _0)`) are set to `true` and used by the assembler later.

A `sparse` assembler is used to arrive at the stiffness matrix and the external load vector using the
finite element requirements as described [here](../01_framework/feRequirements.md#fe-requirements).

```cpp
auto sparseFlatAssembler = SparseFlatAssembler(fes, dirichletValues);
auto req = fe.createRequirement();

auto fextFunction = [&](auto &&lambdaLocal, auto &&dLocal) -> auto & {
  req.insertGlobalSolution( dLocal)
      .insertParameter( lambdaLocal);
  return sparseFlatAssembler.vector(req,VectorAffordance::forces);
};
auto KFunction = [&](auto &&lambdaLocal, auto &&dLocal) -> auto & {
  req.insertGlobalSolution( dLocal)
      .insertParameter( lambdaLocal);
  return sparseFlatAssembler.matrix(req,MatrixAffordance::stiffness);
};
```

The `SparseLU` package from the Eigen library is used to solve the linear system of equations.

For post-processing, the function `Dune::Functions::makeDiscreteGlobalBasisFunction()` is used to create a function for
the displacements and pressure using the basis functions and the nodal values. `Dune::VTKWriter` is used to write
the `*.vtu` files. The results can then be plotted, for example, using [Paraview](https://www.paraview.org/).

```cpp
auto disp
    = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(subspaceBasis(basis.flat(), _0), d);
auto pressure = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis.flat(), _1), d);
Dune::VTKWriter vtkWriter(gridView, Dune::VTK::nonconforming);
vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
vtkWriter.addVertexData(pressure, Dune::VTK::FieldInfo("pressure", Dune::VTK::FieldInfo::Type::scalar, 1));
vtkWriter.write("iks003_incompressibleLinearElasticity");
```

## Takeaways

- `Ikarus::AutoDiffFE` can be used to arrive at the stiffness matrix and external load vector from the energy function.
- Easier implementation of mixed finite elements is possible due to the composite basis feature from Dune.
- Helper functions are included to switch between the Lame parameters.
- Grids from Dune can be directly incorporated within the Ikarus framework.
- Sparse assembler can be used to construct the global stiffness matrices and load vectors.
- Solvers from the Eigen library can be used to solve the linear system of equations.
- Post-processing can be done via Paraview after writing the `*.vtu` files using `Dune::VTKWriter`.
