# Geometry type and description

## Geometry type
{{ inputcpp('src/include/ikarus/Geometries/GeometryType.cpp',False,6,18) }}

The GeometryType represents all possible types of geometry as listed above.
It is used for example in the construction of the grid (see ToDo) and in the
[grind entity interface](theoryGrid.md#interface-of-grid-entity).

## Geometry description
As described in the [interface of the grid entity on the grid theory page](theoryGrid.md#interface-of-grid-entity),
a grid entity has to be able to provide a geometry description, i.e. it has to
return an object which satisfies the geometry interface described below.

### Geometry interface
ToDo: description of the geometry interface

### Geometry implementation
ToDo: description of the available implementations
