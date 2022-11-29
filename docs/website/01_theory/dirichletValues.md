# Dirichlet values
##  Introduction
In finite element problems it is essential to incorporate Dirichlet boundary condition, i.e. fixing part of the solution to some precribed value.
Let's assume we have a problem as

$$
 \boldsymbol{u} = \boldsymbol{g} \quad \text{on } \Gamma_\mathrm{D}.
$$
Here \(  \boldsymbol{u} \) is our solution field and \(\boldsymbol{g}\) is the prescription of it on the Dirichlet boundary \(\Gamma_\mathrm{D}\).

For the discrete algebraic problem this translates to fixing the values of \(u_i\) to some value where \(  \boldsymbol{u}^h = \sum_i N^i u_i \), where
\(  N^i \) is the $i$-th ansatz function and $u_i$ is the solution value at node $i$.

The insertion of several functions $\boldsymbol{g}$ is done in the class `#!cpp Ikarus::DirichletValues`.
##  Interface
The interface of the `#!cpp Ikarus::DirichletValues` is represented by the following code snippet.
```cpp
Ikarus::DirichletValues dirichletValues2(basis); // (1)
void fixBoundaryDOFs(f); // (2)
void fixDOFs(f); // (3)
const auto& basis() const; // (4)
bool isConstrained(std::size_t i) const; // (5)
auto fixedDOFsize() const; // (6)
auto size() const ; // (7) 
```

1. Create class by inserting a [global basis](globalBasis.md)
2. Accepts functor to fix boundary degrees of freedom. "f" is  a functor that will be called with the boolean vector of fixed boundary
 degrees of freedom and the usual arguments of `Dune::Functions::forEachBoundaryDOF`, see Dune book page 388
3. The more general version of `fixBoundaryDOFs`. Here the user should provide a functor that accepts a basis and the correspondign dirichlet degrees of freedom boolean vector
4. Returns the underlying basis
5. Indicates, if a passed degree of freedom is fixed
6. Returns the number of fixed scalar values
7. Returns the number of all dirichlet degrees of freedom

\bibliography