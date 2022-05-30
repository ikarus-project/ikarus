# Global basis

In finite element simulations, we often talk about elements and meshes. What we commonly refer to as "element"
consists of various aspects with different tasks, e.g.  

These elements always connected and disconnected in different ways. 
These connections depend on the underlying basis that is assumed for the solution fields.
The basis functions have usually a local support. E.g. simple 1-D linear Lagrange basis span over two elements.
The connection relation can be encoded in the common vertex node. If we assume higher order 2-D Lagrangian basis function this connection can be associated to a common edge.
There are also ansatz function that only have support within one element. These function are sometimes called bubble-functions.
In the context of discontinuous Galerkin methods the elements are not connected at all.
As last example, if we consider B-Spline basis functions the association of the connection between elements to geometric entities such as edges, vertices fails.

Nevertheless, all this connection information is needed to assemble the global systems matrices. 
However, the global basis needs to provide indices that encode this connectivity depending on the give finite element
This is quite different for different basis. There exists not only one global base.

We are relying here one the basis defined by dune. Thus they use the same interface.
For details see [@sander2020dune] Chapter 10.

\bibliography