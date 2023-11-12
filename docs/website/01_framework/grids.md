# Grids

In the context of the finite element method, the terms "elements" and "meshes" are often used. Each *element* is linked
with certain attributes. The finite element (FE) code should thus be able to evaluate these attributes.
Here are some examples of attributes and the capabilities expected from the FE code:

- *element number* - provide a unique identifier for each element
- *element shape in physical space, shape functions, ...* - provide a description of the geometry
- *element mass matrix, element stiffness matrix, element internal force vector, ...* - provide quantities with physical meaning
- and more ...

In the FE code, there is not one single class that serves to provide all these attributes to an element.
Different goals are achieved by different classes. A group of elements along with their attributes connected to
discretize the actual physical space is called the *mesh*.

In Ikarus and i.e. in Dune, the description of the grid is decoupled from the description of the mesh to preserve its physical meaning.
This helps to have the flexibility of discretizing a simple square-shaped grid with different basis functions like
Langrange or NURBS bases with different polynomial orders. This helps to attain higher levels of abstraction and fewer
iterations of modifying the existing code while studying different grids or the effects of different bases.

For the notions of grids, grid entities, and grid factories, the definitions of Dune are utilized. For details, see
[@sander2020dune] Chapter 5.
All the grids that satisfy the `dune::grid` interface can be used within the Ikarus framework.
This [link](https://www.dune-project.org/doc/grids/) provides an overview of the available `dune::grid` modules.

There also exists an IGA-based grid called [dune-iga](https://github.com/rath3t/dune-iga) to perform
isogeometric analysis (refer [@cottrellIsogeometricAnalysisIntegration2009d]).

It is important to note that the grid only provides geometric information and their relationship to their neighbors.
Even though geometry is typically constructed by some ansatz functions, grids do not provide this information to the
user because some global bases provide ansatz functions for the solution fields. Thus, the user can choose if the
problem should be formulated using the iso-parametric concept, i.e., the same ansatz functions for geometry and
solution fields, or if the ansatz functions for the solution should be independent. The basis defined by Dune is directly
used here. Thus, they use the same interface. For details, see[@sander2020dune] Chapter 10.
\bibliography
