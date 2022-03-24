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
[the implementation of physical quantities can be found here](FiniteElements.md).

On this page, we will go through it using the following example:
![umlDiagram](diagrams/UMLGrid.drawio)

## Grid
The grid is a collection of grid entities. In the example above, the grid consists of three surfaces, 
ten edges and eight vertices, i.e. it consists of 21 grid entities. Since it is some work to construct all
these quantities and their relations, there is the [grid factory](theoryGrid.md#grid-factory) 
which does this job.

### Interface of the grid
- `leafGridView()`: returns a [grid view](theoryGrid.md#grid-view), 
i.e. an object which can iterate over the grid. 

### Available grid implementations
All grids that satisfy the dune::grid interface can be used. For an overview of the available dune::grids, we refer to [link](https://www.dune-project.org/doc/grids/).
Additionally there exists an iga grid [dune-iga](https://github.com/rath3t/dune-iga).

## Grid entity
For the interface of grid entities we refer to [@sander2020dune] Chapter 5.3.

## Grid factory
To construct a grid, a grid factory can be used. To construct a grid, vertices and element
definitions are inserted into the factory. The grid is then constructed by the `createGrid()` function.

### Interface of the grid factory
- `insertVertex(const Dune::FieldVector<double, dimensionworld>&)`: Vertices are inserted as a Dune::FieldVector
  of doubles. Its size is equal to the dimensions of the world (e.g. 2 if it is a 2d simulation etc.)
- `insertElement(Dune::GeometryType type, std::span<size_t> vertices)`: An element is defined
 by its geometrical type and the vertex numbers. For the ordering of the node numbers, see [@sander2020dune] Fig. 5.13.

## Grid view
The interface of a grid view consists of four free functions. Each of them provides a span of
certain grid objects:

- `vertices(GridView)`: returns a span of all vertices in this grid
- `edges(GridView)`: returns a span of all edges in this grid
- `surfaces(GridView)`: returns a span of all surfaces in this grid
- `volumes(GridView)`: returns a span of all volumes in this grid

\bibliography