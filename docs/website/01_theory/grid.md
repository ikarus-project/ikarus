<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de

SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Description of the grid

In finite element simulations, we often talk about elements and meshes. What we commonly refer to as "element"
consists of various aspects with different tasks, e.g.  

- provide a unique identifier (element number)
- provide a description of the reference geometry (element shape in physical space, shape functions, etc.)
- provide quantities with physical meaning (mass matrix, stiffness matrix, internal force vector, etc.)
- ...

In the code, there is not one single class which performs all these tasks. Different tasks are performed by 
different classes, which are described in the following. Especially, the description of the geometry is
decoupled from the task to provide physical meaning. The following content is only about the 
**description of the element geometry**. Details on 
[the implementation of physical quantities can be found here](finiteElements.md).

For the notion of grids, grid entities and grid factories we rely on the definitions of dune. For details, see
[@sander2020dune] Chapter 5.

### Available grid implementations
All grids that satisfy the dune::grid interface can be used. For an overview of the available dune::grids, we refer to [link](https://www.dune-project.org/doc/grids/).
Additionally, there exists an iga grid [dune-iga](https://github.com/rath3t/dune-iga).

\bibliography