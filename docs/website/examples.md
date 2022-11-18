# Examples

In order to understand several features of Ikarus, a set of examples are provided within the finite element framework.
These can be found in [IkarusExamples](https://github.com/IkarusRepo/IkarusExamples). The installation and execution 
methodologies are briefly commented in the [README](https://github.com/IkarusRepo/IkarusExamples/blob/main/README.md) file 
of the repository. Each example is given a unique identification in the beginning of the file name of the form `iksXXX`.
This unique identification is also used in the following instead of the complete `*.cpp` file name. The auxiliary files 
to the examples like `*.msh`, `*.geo` or `*.parset` are available in `../../src/testFiles/`. 
In order to add a new example, create a pull request with your executable file in the repository IkarusExamples and in 
parallel update the documentation here, see [How to contribute](https://ikarusrepo.github.io/codeStyle/) and 
[How to edit](https://ikarusrepo.github.io/documentation/howToEdit/) for more information.  
The available examples are described in the following.

## Cantilever beam with point load
The example `iks001_cantileverBeam_oneDGrid.cpp` shows a simple implementation of a one dimensional Timoshenko beam which is clamped on the left 
hand side. A point load is applied on the right hand side of the beam. It uses `Dune::OneDGrid` to generate the required 
grid. A simple implementation is shown here where the stiffness matrices are assembled explicitly. Advanced 
implementations of matrix assembly and other features of Ikarus is showcased in the other examples.

## Compute the value of $\pi$
The examples `iks002_compute-pi.cpp` and `iks003_compute-pi.cpp`shows the calculation of $\pi$ by computing the area 
and circumference of a unit circle. These examples help to understand the `Grid` module from Dune and the refinement techniques it 
brings. The example `iks002` shows that a global refinement doesn't refine the number of grid entities on the boundary 
of the circle, which leads to a poor approximation of $\pi$ while comparing with the circumference of the circle. 
On the other hand, `iks003` shows how elements on the boundaries can be marked and refined, thereby resulting in an 
accurate approximation of $\pi$.

## Compression of an incompressible rubber block
`iks004_incompressible_LinearElasticity.cpp` uses a finite element technology with displacement and pressure as 
independent degrees of freedom to simulate the compression of an incompressible rubber block. The potential energy for such a system is defined in the 
`calculateScalarImpl(const FERequirementType &par, const Eigen::VectorX<ScalarType> &dx)` function in `struct` 
named `Solid`. This function uses the principles of automatic differentiation to provide the stiffness matrices and 
other necessary quantities to provide a static structural analysis.   

## Plate subjected to a surface load
Kirchhoff type plate element is implemented in `iks005_kirchhoff-plate.cpp` using the automatic differentiation 
technique as commented before. The basis used for discretization is a NURBS basis from the `dune-iga` module.
The problem is solved and convergence plots are created by comparing the solutions to available analytical solutions for 
simply supported and clamped boundaries.

## Newton-Raphson method
`iks006_newtonRaphson.cpp` shows a basic example of the Newton-Raphson method to solve a non-linear set of equations. 
A function which shows the algorithm explicitly is provided and another function which is implemented in Ikarus is 
demonstrated. The function which depicts the Ikarus implementation uses a 
[non-linear operator](https://ikarusrepo.github.io/01_theory/nonlinearOperator/) to 
perform the Newton-Raphson iterations. A logger can also be subscribed to in order to observe the residual norms, 
for instance.

## Non-linear for 2D solids
Again, automatic differentiation based implementation is used to perform a non-linear analysis for a 2D block in 
`iks007_nonlinear2Dsolid.cpp`. Various methods to obtain a 2D grid via Dune is also shown in the commented section in 
the beginning. Python is used to provide a Neumann boundary condition providing a demonstration for the usage of a 
Python-based code within the Ikarus framework. Load control method is chosen as the desired control routine and 
Newton-Raphson (or Trust region methods) are used to solve the non-linear problem itself.

## Von-Mises stress calculation for truss systems
`iks008_vonmises_truss.cpp` shows a way to use the tools and features mentioned in the previous examples to calculate 
and post-process the Von-Mises stresses in truss systems.

## Cook's membrane
The Cook's membrane problem adapted from the paper[@cook_improved_1974] is implemented in examples
`iks009_cook_membrane.cpp` and `iks010_cook_membrane_convergence.cpp`. This problem can be solved not only with
structured meshes provided, but also with unstructured and triangular meshes. The input parameters like material and grid 
parameters are read from the file `cook.parset`. The problem can be solved also with the standard planar solid element, 
or with enhanced assumed strain elements. For more details on the element technologies, refer the 
[documentation](https://ikarusrepo.github.io/01_theory/finiteElements/). `iks009` solves the problem for a chosen 
finite element type whereas `iks010` solves the problem with a set of existing finite elements and compares the 
convergence rates. 