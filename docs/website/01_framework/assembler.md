# Assembler

The purpose of an assembler is to assemble local quantities (local stiffness matrix, local force vector, local energy, etc.)
by looping over finite elements and thereby arriving at a global structure. This page describes
the available assemblers and how they can be used.

Each of the assemblers is constructed as follows:

```cpp
AssemblerName(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichletFlags)
```

- `basis` is the basis that was used to construct the finite elements. The implementation of the bases involves four different
  strategies: `BlockedLexicographic`, `BlockedInterleaved`, `FlatLexicographic` and `FlatInterleaved`. These strategies are to be
  considered while creating the assembler. For further information on these strategies, refer to Chapter 10 of DUNE[@sander2020dune].
- `fes` is a container that contains all the finite elements that should be assembled.
- `dirichletFlags` is a `#!cpp std::vector<bool>` type. The `i`-th degree of freedom is fixed when `#!cpp dirichletFlags[i] = true`.
  When a reduced matrix or vector is chosen, the corresponding row and column entries are removed.

## Base class for flat assemblers

The FlatAssemblerBase is the base class for all assemblers currently available. All other assemblers inherit from this one,
i.e., their interface includes the following functions:

```cpp
size_t size() // (1)!
size_t reducedSize() // (2)!
auto &finiteElements() const // (3)!
Eigen::VectorXd createFullVector(const Eigen::VectorXd &reducedVector) // (4)!
size_t constraintsBelow(size_t i) // (5)!
bool isConstrained(size_t i) // (6)!
size_t estimateOfConnectivity() // (7)!
```

1. Returns the number of degrees of freedom.
2. Returns the number of degrees of freedom that are not constrained by a Dirichlet boundary condition.
3. Returns a reference to the finite element container, which was passed to the assembler.
4. Gets a reduced vector and returns a full vector. Entries corresponding to fixed dofs are set to 0. The values of the other entries are
    obtained from the reduced vector.
5. Indicates how many degrees of freedom {0,1,...i-1} are fixed.
6. Indicates whether the degree of freedom `i` is fixed.
7. Returns 8x the number of grid elements, which is an estimate for the connectivity. It can be used to allocate vectors.

## Scalar assembler

It has the capabilities of [FlatAssemblerBase](#flatassemblerbase) plus one additional function:

```cpp
double& getScalar(const RequirementType& fErequirements)
```

This assembler can be used when only a scalar quantity is of interest and the assembly of matrices or vectors is irrelevant.
The available requirements are explained on the [FE requirements page](feRequirements.md).
`dirichletFlags` is not used in this assembler.

It assembles the requested scalar quantity. A call to this function could look like this:

```cpp
ScalarAssembler myAssembler(...) // (1)!
const auto& K = myAssembler.getScalar(feRequirements) // (2)!
```

1. Represents the construction of the desired assembler.
2. To learn more about the available alternatives to `energy` and how this works, read the [FE requirements](feRequirements.md) page.

## Flat vector assembler

It has all the features of [ScalarAssembler](#scalarassembler) plus more, like:

```cpp
Eigen::VectorXd& getRawVector(const FERequirementType &feRequirements)
Eigen::VectorXd& getVector(const FERequirementType& feRequirements)
Eigen::VectorXd& getReducedVector(const FERequirementType& feRequirements)
```

The first one returns a vector without considering the boundary conditions.
The remaining, as the name suggests, returns a full vector or a reduced vector considering boundary conditions.
They work in the same way as the scalar assembly functions of [ScalarAssembler](#scalarassembler).
The available FE requirements are explained on the [FE requirements](feRequirements.md) page.

## Flat sparse assembler

It offers the functions of [VectorFlatAssembler](#vectorflatassembler) plus more, like:

```cpp
Eigen::SparseMatrix<double> &getRawMatrix(const FERequirementType &feRequirements)
Eigen::SparseMatrix<double> &getMatrix(const FERequirementType &feRequirements)
Eigen::SparseMatrix<double> &getReducedMatrix(const FERequirementType &feRequirements)
```

A sparse matrix is returned.
They work in the same way as the vector assembly functions of [VectorFlatAssembler](#vectorflatassembler).
The available FE requirements are explained on the [FE requirements](feRequirements.md) page.

## Flat dense assembler

The only difference between the [SparseFlatAssembler](#sparseflatassembler) and the DenseFlatAssembler is that the
DenseFlatAssembler returns a dense matrix.

```cpp
Eigen::MatrixXd &getRawMatrix(const FERequirementType &feRequirements)
Eigen::MatrixXd &getMatrix(const FERequirementType &feRequirements)
Eigen::MatrixXd &getReducedMatrix(const FERequirementType &feRequirements)
```
