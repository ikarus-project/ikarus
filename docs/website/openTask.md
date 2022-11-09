# Open tasks
Thanks for your interest in contributing to this code base.
If your are interested the following task are vacant.


### Local functions
* Implementing a unit normal field function[^KLNote] and its derivatives w.r.t. its coefficients \( \boldsymbol{x}_i \)

    \[ 
    \boldsymbol{n} = \frac{\boldsymbol{a}_1 \times \boldsymbol{a}_2}{||\boldsymbol{a}_1 \times \boldsymbol{a}_2||}, \quad \text{with } \boldsymbol{a}_{\alpha} = \sum_{i=1}^n N^i_{,\alpha}(\boldsymbol{\xi}) \boldsymbol{x}_i
    \] 

    To implement these see [link](01_theory/localFunctions.md#how-to-implement-your-own-local-functions).

* Support second derivatives
* Add \( \operatorname{div} \) and \( \operatorname{curl} \) wrapper

### Controlroutines
* Dynamics (Explicit/ implicit time stepping)

### Controlroutines addons
* Extended systems
* Inhomogeneous dirichlet boundary conditions wrapper class

### Finite element helper
* Implement default implemented mass matrix

### Finite elements
* Nonlinear Reissner-Mindlin shell [@muller2022consistent]
* Kirchhoff-Love shell
* 3D-Beam
* Implement forces and stiffness matrix of `NonLinearElasticityFE`

### Local Basis 
* Support second derivatives


### Addons
* Add Python binding [pybind11](https://github.com/pybind/pybind11)
* Add [Muesli](https://materials.imdea.org/muesli/)

[^KLNote]: This is usually needed for a Kirchhoff-Love shell implementation, see [@kiendlKLshell].


!!! note  "Code style"
For details on our code style we refer to [Link](codeStyle.md).


\bibliography 