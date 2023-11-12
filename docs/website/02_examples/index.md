# Examples

In order to understand several features of Ikarus, a set of examples is provided within the finite element framework.
These can be found at [IkarusExamples](https://github.com/ikarus-project/ikarus-examples). The installation and execution
methodologies are briefly discussed in the [README](https://github.com/ikarus-project/ikarus-examples/blob/main/README.md) file
of the repository. Each example is given a unique identification in the beginning of the file name of the form `iksXXX`.
This unique identification is also used in the following instead of the complete `*.cpp` file name. Auxiliary files
for the examples, such as `*.msh`, `*.geo`, or `*.parset` can be found in `../../src/testfiles/`.
In order to add a new example, create a pull request with your executable file in the repository IkarusExamples and, in
parallel, update the documentation here. See [How to Contribute](../03_contribution/codeStyle.md) and
[How to Edit](../03_contribution/howToEdit.md) for more information.

The available examples are:

 | Identification | Name of the example                                                           |
|:---------------|:------------------------------------------------------------------------------|
| iks001         | [Compute the value of $\pi$](computePi.md)                                    |
 | iks002         | [Cantilever beam with point load](cantileverBeam.md)                          |
| iks003         | [Deformation of an incompressible rubber block](incompressibleRubberBlock.md) |
| iks004         | [Plate subjected to a surface load](kirchhoffPlate.md)                        |
| iks005         | [Newton-Raphson method](newtonRaphsonMethod.md)                               |
| iks006         | [Non-linear Elasticity for 2D solids](nonLinearElasticity.md)                 |
| iks007         | [Von-Mises truss](vonMisesTruss.md)                                           |
| iks008         | [Cook's membrane](cooksMembrane.md)                                           |
