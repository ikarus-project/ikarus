# Grid tutorials
This examples explains how to use and construct grids.
The syntax and functionality is closely related to grid objects from the [DUNE Numerics project](https://dune-project.org).  


### SimpleGrid
This tutorial explains how to use `SimpleGrid<dim,wdim>`. 

`SimpleGrid<dim,wdim>` uses two mandatory non-type template integers. The first integer denotes 
the largest possible dimension of a `GridEntity` and the second integer denotes the embedding 
space of the grid. E.g. `SimpleGrid<2,3>` contains surfaces,edges,vertices living in three dimensions.

We go step by step through the following code. The code can be found in `examples/src/drawgrid.cpp`

{{ inputcpp('examples/src/drawgrid.cpp',True) }}

First we choose, which grid type we want to create. In this case we create a grid of 3d elements
living in 3d space. Therefore, we create a factory for this type.
{{ inputcpp('examples/src/drawgrid.cpp',False,7,9) }}
After that we fill a vector with nine vertex positions. The type of the vector of coordinates is the three dimensional
`double`vector `Eigen::Vector3d` of the [Eigen library](https://eigen.tuxfamily.org/).
{{ inputcpp('examples/src/drawgrid.cpp',False,9,22) }}
In the following we use the first functionality of the gridfactory which is the `insertVertex()` 
function which allows to insert vertices into grid. We insert all vertices previously defined.
{{ inputcpp('examples/src/drawgrid.cpp',False,21,24) }}
After that, we define elements. The first elements is a linear hexahedron defined by vertices 1-7:
{{ inputcpp('examples/src/drawgrid.cpp',False,24,29) }}
Furthermore, we define a linear tetrahedron:
{{ inputcpp('examples/src/drawgrid.cpp',False,29,32) }}
Finally we create the grid and a grid view to access the grid:
{{ inputcpp('examples/src/drawgrid.cpp',False,32,36) }}
To visualize the grid, we use a draw function:
{{ inputcpp('examples/src/drawgrid.cpp',False,37,38) }}
<figure>
  <img src="/images/simpleGrid.svg" width="600" />
  <figcaption>Grid containing a hexahedron and a tetrahedron. In red the vertices are shown. The grid shows a non-conforming mesh.</figcaption>
</figure>
