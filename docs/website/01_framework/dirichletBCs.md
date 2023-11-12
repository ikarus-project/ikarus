# Dirichlet boundary conditions

## Introduction

In finite element problems, it is essential to incorporate the Dirichlet boundary conditions, i.e.,
prescribing a part of the solution to a fixed value.
Let us consider the following:

$$
 \boldsymbol{u} = \boldsymbol{g} \quad \text{on  } \Gamma_\mathrm{D}.
$$
Here \(  \boldsymbol{u} \) is the solution field with a prescribed value \(\boldsymbol{g}\) on the boundary \(\Gamma_\mathrm{D}\).

For the discrete algebraic problem, this translates to fixing the values of \(u_i\) to \(g_i\) in the approximation
\(  \boldsymbol{u}^h = \sum_i N^i u_i \), where \(  N^i \) is the $i$-th ansatz function.

The handling of such a function $\boldsymbol{g}$ is done by the class `#!cpp Ikarus::DirichletValues`.

## Interface

The interface of `#!cpp Ikarus::DirichletValues` is represented by the following code snippet:

```cpp
Ikarus::DirichletValues dirichletValues2(basis); // (1)!
void fixBoundaryDOFs(f); // (2)!
void fixDOFs(f); // (3)!
const auto& basis() const; // (4)!
bool isConstrained(std::size_t i) const; // (5)!
auto fixedDOFsize() const; // (6)!
auto size() const ; // (7)!
```

1. Create class by inserting a global basis, [@sander2020dune] Chapter 10.
2. Accepts a functor to fix boundary degrees of freedom. `f` is  a functor that will be called with the boolean vector of fixed boundary.
 degrees of freedom and the usual arguments of `Dune::Functions::forEachBoundaryDOF`,  as defined on page 388 of the Dune
   [@sander2020dune] book.
3. A more general version of `fixBoundaryDOFs`. Here, a functor is to be provided that accepts a basis and the corresponding boolean
   vector considering the Dirichlet degrees of freedom.
4. Returns the underlying basis.
5. Indicates whether the degree of freedom `i` is fixed.
6. Returns the number of fixed degrees of freedom.
7. Returns the number of all dirichlet degrees of freedom.

\bibliography
