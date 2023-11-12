# Compute the value of $\pi$

## Description

The example `iks001_computePi.cpp` shows the calculation of $\pi$ by computing the area
and circumference of a unit circle. This example helps to understand the `Grid` module from Dune and the refinement techniques it
brings. The example shows that a global refinement doesn't refine the number of grid entities on the boundary
of the circle, which leads to a poor approximation of $\pi$ when comparing it with the area and circumference of the circle.
On the other hand, it also shows how elements on the boundaries can be marked and refined, thereby resulting in an
accurate approximation of $\pi$.

## Code highlights

This example contains the following two functions:

```cpp
boundaryUnawareRefinedCircle();
boundaryAwareRefinedCircle();
```

demonstrating the usage of certain attributes of the grid module from Dune[@sander2020dune].
Firstly, we explain the function `#!cpp boundaryUnawareRefinedCircle();`.

Here, a `circleCoarse.msh` file, which was created using [Gmsh](https://gmsh.info/), is read using the `Dune::GmshReader` and a
`Dune::ALUGrid` object is created as shown below:

```cpp
constexpr int gridDim = 2;  // (1)
using Grid            = Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>;
auto grid             = Dune::GmshReader<Grid>::read("auxiliaryFiles/circleCoarse.msh", false);
auto gridView         = grid->leafGridView();  // (2)
```

It is to note that dune-grid only supports Version 2 of the Gmsh format. The `#!cpp draw(gridView);` functionality is
included within the Ikarus framework to quickly draw grids and verify them for any major errors.
The function `#!cpp grid->globalRefine(1);` is invoked to refine the grid on a global level. This means that, if we have
a single square-shaped 4-node quadrilateral element, the `globalRefine(1)` function will bisect the element once in
either direction and thereby result in a grid with 2 elements in either direction. This type of refinement is done for triangular elements in
the following, and the area of the circle is compared to the value of $\pi$. The area of the circle itself is obtained by
summing the volumes (area in 2D terms) of individual elements.

```cpp
double area = 0.0;
for (int i = 0; i < 3; ++i) {
  area = 0.0;
  grid->globalRefine(1);
  auto gridViewRefined = grid->leafGridView();
  std::cout << "This gridview contains: ";
  std::cout << gridViewRefined.size(0) << " elements" << std::endl;
  for (auto &element : elements(gridViewRefined)) {
    area += element.geometry().volume();
  }
  std::cout << std::setprecision(10) << "Area: " << area << " Pi: " << std::numbers::pi << std::endl;
}
```

The area of individual elements is also written in a `*.vtu` file using the module `Dune::VTKWriter`. More explanations
of this are found in subsequent examples. In this example, we restrict ourselves to the grid refinement strategies and
accessing grid entities from Dune.

It is also possible to check if the element has any edges at the boundary using the method `#!cpp element.hasBoundaryIntersections()`.
If this is true, the edges intersecting with the boundaries can be extracted, and the `volume` method on this `intersection`
object would then return the length of the edge in the 2D case. The circumference of this unit circle is then computed as shown below:

```cpp
double circumference = 0.0;
for (auto &element : elements(gridView))
  if (element.hasBoundaryIntersections())
    for (auto &intersection : intersections(gridView, element))
      if (intersection.boundary()) circumference += intersection.geometry().volume();
```

It is observed that even though the grid is globally refined in every loop, the computed value of $\pi$ never converges.
This is due to the fact that the triangular grid doesn't change its area with refinement. Even though the grid has many elements with every refinement,
it doesn't contain enough information about the boundary of the circle to predict an accurate value of $\pi$. In order to refine
the boundary entities correctly, the `#!cpp boundaryAwareRefinedCircle()` function is formulated and explained next.

In this function, the grid is created explicitly and isn't read from an external `*.msh` file. The corners are
calculated manually, which is followed by the insertion of vertices, elements, and boundary segments as shown below:

```cpp
Dune::GridFactory<Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>> gridFactory;
Eigen::Vector2d v(1, 0);
std::array<FieldVector<double, 2>, 6> corners0;
Eigen::Rotation2D<double> R;
R.angle() = 0.0;
for (auto &corner : corners0) {
  Eigen::Vector2d a = R * v;
  corner[0]         = a[0];
  corner[1]         = a[1];
  R.angle() += 60.0 / 180.0 * std::numbers::pi;
}

gridFactory.insertVertex({0, 0});
gridFactory.insertVertex(corners0[0]);
gridFactory.insertVertex(corners0[1]);
gridFactory.insertVertex(corners0[2]);
gridFactory.insertVertex(corners0[3]);
gridFactory.insertVertex(corners0[4]);
gridFactory.insertVertex(corners0[5]);

gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 1, 2});
gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 2, 3});
gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 3, 4});
gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 4, 5});
gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 5, 6});
gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 6, 1});

/// Create boundary segments which map the boundaries onto the unit circle
gridFactory.insertBoundarySegment({1, 2}, std::make_shared<UnitCircleBoundary>(corners0[0], corners0[1]));
gridFactory.insertBoundarySegment({2, 3}, std::make_shared<UnitCircleBoundary>(corners0[1], corners0[2]));
gridFactory.insertBoundarySegment({3, 4}, std::make_shared<UnitCircleBoundary>(corners0[2], corners0[3]));
gridFactory.insertBoundarySegment({4, 5}, std::make_shared<UnitCircleBoundary>(corners0[3], corners0[4]));
gridFactory.insertBoundarySegment({5, 6}, std::make_shared<UnitCircleBoundary>(corners0[4], corners0[5]));
gridFactory.insertBoundarySegment({6, 1}, std::make_shared<UnitCircleBoundary>(corners0[5], corners0[0]));

auto grid     = gridFactory.createGrid();
auto gridView = grid->leafGridView();
```

To refine the grid entities living in the boundary, the elements that have intersections with the boundary are
first marked and then refined, as shown in the code below:

```cpp
for (const auto &ele : elements(grid->leafGridView())) {
  if (ele.hasBoundaryIntersections()) grid->mark(1, ele);
}
grid->preAdapt();
grid->adapt();
grid->postAdapt();
auto gridViewRefined = grid->leafGridView();
```

Now, the calculation of area and circumference to determine the value of $\pi$ converges correctly.

## Takeaways

- `#!cpp Dune::GmshReader` can be used to read `*.msh` files to import grids.
- Grid entities can be marked and locally refined, or the grid can be globally refined.
- Grids can also be explicitly created by inserting vertices, elements, and boundary segments.
- `#!cpp element.hasBoundaryIntersections()` can be used to check if an element has any intersections with the boundaries.
