# Grid tutorials
This examples explains how to use and construct grids.
The syntax and functionality is closely related to grid objects from the [DUNE Numerics project](https://dune-project.org).  


### SimpleGrid
There, exists a simple grid implementation called `SimpleGrid<dim,wdim>`. It can be used for up to three world space 
dimension and it provides an unstructured grid with arbitrary element types of the same dimension.
It uses two mandatory non-type template integers. The first integer denotes the largest possible dimension of a `GridEntity` 
and the second integer denotes the embedding space of the grid. E.g. `SimpleGrid<2,3>` contains surfaces,edges,vertices living in three dimensions.

![umlDiagram](../diagrams/UMLGrid.drawio)
### SimpleGrid example
We go step by step through the following code.

{{ inputcpp('examples/src/drawgrid.cpp',True) }}

First we choose, which grid type we want to create and fill with entities using the type alias `Grid`.
{{ inputcpp('examples/src/drawgrid.cpp',False,7,9) }}
 After this we create a factory for this type.
{{ inputcpp('examples/src/drawgrid.cpp',False,9,10) }}
After that we fill a vector with nine vertex positions. The type of the vector of coordinates is the three dimensional
`double`vector `Eigen::Vector3d` of the [Eigen library](https://eigen.tuxfamily.org/).
{{ inputcpp('examples/src/drawgrid.cpp',False,10,22) }}
In the following we use the first functionality of the gridfactory which is the `insertVertex()` function which allows to insert vertices into grid
{{ inputcpp('examples/src/drawgrid.cpp',False,22,24) }}
<figure>
  <img src="/images/simpleGrid.svg" width="600" />
  <figcaption>Grid containing a hexahedron and a tetrahedron. In red the vertices are shown. The grid shows a non-conforming mesh.</figcaption>
</figure>

## Interfaces
### Grid factory
### Grid
The grid itself does not have a particular interface, except that it has to provide function to return a grid view.

| Grid Entity Interface        ||
| :------------ | :-----------: |
| `#!cpp GridViewType leafGridView()`     |
| `#!cpp GridViewType levelGridView(int levl)`     |

### Grid entities
To stay compatible with `Dune::GridEntities` we provide a subset of the interface for the grid entities. 
Furthermore, we demand that they need to provide a unique identitfier.


```cpp
                   // Return a unique id of this entity
size_t             getID();   
                   // The refinement level to which the entity belongs
int                level();   
                   // Geometric realization of this entity, see
GeometryType       geometry();   
                   // The type of the entity, e.g. line, vertex, quadrilateral,...
Dune::GeometryType type();
                   // Number of subentities, e.g. a line has two vertices
unsigned int       subEntities();
```
For the interface of `GeometryType` see [Geometry](../theory/theoryGrid.md).


!!! warning
    It is tempting to attach any physical meaning to the grid entities. Nevertheless, the only purpose of the grid and grid entities
    is to provide the connectivity between different grid entities and an unique identifier to construct later on degrees of freedom indices from it.


## Formula test

$$
\operatorname{ker} f=\{g\in G:f(g)=e_{H}\}{\mbox{.}}
$$