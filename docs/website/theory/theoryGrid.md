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
[the implementation of physical quantities can be found here](theoryFiniteElement.md).

The grid description will be explained using the following example
![img.png](images/sampleProblemGrid.png)

## Grid
The grid is a collection of grid entities. In the example above, the grid consists of three surfaces, 
ten edges and eight vertices, i.e. it consists of 21 grid entities. Since it is some work to construct all
these quantities and their relations, there is the [grid factory](theoryGrid.md#grid-factory) 
which does this job.

### Grid interface
The grid only has to implement one function: `leafGridView()` which returns a [grid view](theoryGrid.md#grid-view), 
i.e. an object which can iterate over the grid. 

### Grid implementation
There is currently one implementation which is `SimpleGrid`. ToDo: describe implementation details

## Grid entity
A grid entity provides all information related to the element geometry. 

### Interface of grid entity

A grid entity has the following properties:
- It has a certain
  [geometrical type (see below)](theoryGrid.md#geometry-type), e.g. it can be a vertex, a linear line defined by
  two point or a quadrilateral with linear edges. The type can be obtained by calling `type()`
- It has sub-entities of lower dimension. In the example above, surface S1 has the following sub-entities:
  - four edges (1,2,4,6)
  - four vertices (1,2,4,5)
  
  The number of each of these sub-entities can be obtained from `subEntities(codim)` where codim means ???
- The subEntities itself can be accessed by ???  
- It can provide a geometry description. By calling the `geometry()`function, an object is return which
  the geometry interface. Further details about this interface and what the returned object is able to do
  can be found on the [geometry theory page](theoryGeometry.md).
- It has an unique identifier, which can be obtained by calling the member funciton `getID()`
- It knows its dimension (e.g. a line is 1D, a quadrilateral is 2D) which can be obtained from `mydimension`
- It knows the dimension of the world it lives in, which can be
  obtained from `dimensionworld`. If it is a 1d truss element, it can tell whether it lives in a 1d, 2d or 3d space.
- It knows its `codimension`, which is ???
- ToDo after refactoring this interface, complete the documentation. ToDo unify with UML diagram.

```cpp
  concept GridEntity = requires(GridEntityType gEntity, unsigned int codim) {
    { gEntity.level() } -> std::same_as<int>;
    gEntity.geometry();
    { gEntity.type() } -> std::same_as<Ikarus::GeometryType>;
    { gEntity.subEntities(codim) } -> std::same_as<unsigned int>;
    GridEntityType::codimension;
    GridEntityType::dimensionworld;
    GridEntityType::mydimension;
    typename GridEntityType::Geometry;
  };
```

### Implementation of grid entity
There is currently one implementation of the GridEntity interface available, which is `DefaultGridEntity`. It is
supposed to be used together with the grid implementation `SimpleGrid`.
```cpp
  template <int griddim, int cogriddim, int wdim>
  class DefaultGridEntity {
    // ...
  };
```
It is based on three template parameters `griddim`(dimension of the grid), `cogriddim`(???) and 
`wdim`(dimension of the world).

ToDo: How are the private quantities (e.g. entitiesFathers) constructed?

## Grid factory
To construct a grid, a grid factory can be used. 

### Interface of the grid factory
ToDo: Discussion: grid factory is very specific for the grid which is constructed. Is an interface required?

### Implementation of grid factory
There is currently one implementation of a grid factory which is `SimpleGridFactory`. 
It constructs a `SimpleGrid`. 


## Grid view
TODO
### Interface of grid view
TODO
### Implementation of grid view
TODO