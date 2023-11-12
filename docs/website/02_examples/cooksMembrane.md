# Cook's membrane

## Description

The example `iks008_cooksMembrane.cpp` implements the Cook's membrane problem adapted from the paper[@cook_improved_1974].
This problem can be solved not only with the structured meshes provided, but also with unstructured and triangular meshes.
The input parameters, like material and grid
parameters, are read from the file `iks008_cooksMembrane.parset`. The problem can also be solved with the standard 2D planar solid element
or with enhanced assumed strain elements. For more details on the element technologies, refer to the
[documentation](../01_framework/finiteElements.md). `iks008` solves the problem with a set of existing finite elements and compares the
convergence rates.

## Code highlights

In this example, the standard Q1 finite element and the enhanced assumed strain elements Q1E4, Q1E5, and Q1E7 are used to solve the
Cook's membrane problem.
Convergence studies are done for the vertical displacement in the top right corner, and the assembly time for the stiffness matrix is
also compared.
The [date and time utilities](https://en.cppreference.com/w/cpp/chrono) from the standard C++ library are used to determine the
computation times in this example.
This example not only loops over the different refinement levels for a convergence plot but also loops over the different element types.
In order to avoid re-building the complete code with modifications to certain input parameters, the `parametertreeparser.hh` can be used
from the [dune-common](https://www.dune-project.org/modules/dune-common/) module.
The `Dune::ParameterTree` is used to read the Young's modulus (`E`), the Poisson's ratio (`nu`), and the level of refinement
(`refinement_level`), as shown below:

```cpp
Dune::ParameterTree parameterSet;
Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

const Dune::ParameterTree &gridParameters     = parameterSet.sub("GridParameters");
const Dune::ParameterTree &controlParameters  = parameterSet.sub("ControlParameters");
const Dune::ParameterTree &materialParameters = parameterSet.sub("MaterialParameters");

const double E             = materialParameters.get<double>("E");
const double nu            = materialParameters.get<double>("nu");
const int refinement_level = gridParameters.get<int>("refinement");
```

`argv[1]` is the variable in the argument vector that contains `iks008_cooksMembrane.parset` to read the input parameters.
The file `cook.msh` contains the Cook's membrane problem with a structured grid, whereas `cook_tri.msh` and `cook_unstructured.msh`
contains the same problem with triangular elements and with an unstructured mesh, respectively. The mesh file is read using
`Dune::GmshReader` and `Dune::UGGrid` is used to get the `grid` object.
`easSet` is an `Eigen::Vector` that contains the number of EAS parameters for the four element types.
It is important to note that if the number of EAS parameters is set to zero, the standard Q1 formulation is used.

```cpp
Eigen::Vector<int, 4> easSet;
easSet << 0, 4, 5, 7;
```

The EAS elements are created then, as shown below:

```cpp
auto numberOfEASParameters = easSet(nep); // (1)!
std::vector<Ikarus::EnhancedAssumedStrains<Ikarus::LinearElastic<decltype(basis)>>> fes;
for (auto &element : elements(gridView)) {
  fes.emplace_back(basis, element, E, nu, &volumeLoad, &neumannBoundary, &neumannBoundaryLoad);
  fes.back().setEASType(numberOfEASParameters);
}
```

1. `nep` is the index of the `for`-loop which runs from 0 to 4 here.

The Dirichlet boundary conditions are defined for the left edge, and the Neumann boundary condition on the right edge is defined by the
usage of the dune-python interface.
A sparse assembler is used, and the linear system of equations is solved. The vertical displacement in the top right corner is computed
as shown below and is later stored in a vector:

```cpp
auto dispGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(*basis, D_Glob);
auto localView      = basis.flat().localView();
auto localw         = localFunction(dispGlobalFunc);
double uy_fe        = 0.0;
Eigen::Vector2d req_pos;
req_pos << 48.0, 60.0;
for (auto &ele : elements(gridView)) {
  localView.bind(ele);
  localw.bind(ele);
  const auto geo = localView.element().geometry();
  for (size_t i = 0; i < 4; ++i) {
    if (Dune::FloatCmp::eq(geo.corner(i)[0], req_pos[0]) and Dune::FloatCmp::eq(geo.corner(i)[1], req_pos[1])) {
      const auto local_pos = geo.local(toDune(req_pos));
      uy_fe                = toEigen(localw(local_pos)).eval()[1];
    }
  }
}
```

The datasets are then stored and plotted using [Matplot++](https://github.com/alandefreitas/matplotplusplus). The deformed configuration
is also written using the `Dune::VTKWriter` and can be visualized using Paraview.
Several log information is also displayed in this example using [spdlog](https://github.com/gabime/spdlog).

## Takeaways

- `Dune::ParameterTree` can be used to read input parameters from an external `*.parset` file.
- A linear elastic element of arbitrary dimension can be used to solve the underlying problem. In 2D and 3D, this element can be
  enriched with enhanced assumed strains.
- `spdlog/spdlog.h` can be used to display log information.
- `chrono` library can be used to determine the computation time.
