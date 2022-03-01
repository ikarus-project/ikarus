# Grid tutorials
This examples explains how to use and construct grids.
The syntax and functionality is closely related to grid objects from the [DUNE Numerics project](https://dune-project.org).  


### SimpleGrid
This tutorial explains how to use `SimpleGrid<dim,wdim>`. 

`SimpleGrid<dim,wdim>` uses two mandatory non-type template integers. The first integer denotes 
the largest possible dimension of a `GridEntity` and the second integer denotes the embedding 
space of the grid. E.g. `SimpleGrid<2,3>` contains surfaces,edges,vertices living in three dimensions.

