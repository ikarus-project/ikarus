# Assembler

The purpose of an assembler is to assemble global quantities (global stiffness matrix, global force vector, global energy, ...) 
by looping over finite elements and composing the corresponding local structures to a global structure. This page describes
the available assemblers and how they can be used.

Each of the assemblers is constructed as follows:
```cpp
AssemblerName(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
```

- `basis` is the basis that was used to construct the finite elements. ToDo add comment about FLAT.
- `fes` is a container that contains all the finite elements that should be assembled
- `dirichFlags` is of type `#!cpp std::vector<bool>`. `#!cpp dirichFlags[i] = true` means that degree of freedom i is fixed. 
    The corresponding row / column / entry will be eliminated when you ask for reduced matrix / vector. 

## FlatAssemblerBase
The FlatAssemblerBase is the basis for all assemblers currently available. All other Assemblers inherit from this assembler, 
i.e. they have the functions listed below as well:
```cpp
size_t size() // (1)
size_t reducedSize() // (2)
auto &finiteElements() const // (3)
Eigen::VectorXd createFullVector(const Eigen::VectorXd &reducedVector) // (4)
size_t constraintsBelow(size_t i) // (5)
bool isConstrained(size_t i) // (6)
size_t estimateOfConnectivity() // (7)
```

1. Returns the number of degrees of freedom.
2. Returns the number of degrees of freeedom, which are not constrained by a dirichlet boundary condition.
3. Returns a reference to the finite element container that you gave to the assembler when constructing it.
4. Gets a reduced vector and returns a full vector. Entries corresponding to fixed dofs are set to 0. Values of the other entries are
    obtained from the reduced vector.
5. Tells you how many of the degrees of freedom {0,1,...i-1} are fixed.
6. Tells you if degree of freedom i is fixed
7. Returns 8x the number of grid elements, which is an estimate for the connectivity. It can be used to allocate vectors.


## ScalarAssembler
It has the capabilities of [FlatAssemblerBase](#flatassemblerbase) plus one additional function:
```cpp
double& getScalar(const RequirementType& fErequirements)
```
This assembler can be used when you are only interested in a scalar quantity
and assembling of matrices or vectors is not relevant for you.
The available requirements are explained on the [FE requirements page](feRequirements.md).
`dirichletFlags` is not used in this assembler.

It assembles the requested scalar quantity. A call to this function could look as follows:
```cpp
ScalarAssembler myAssembler(...) // (1)
// other code
const auto& K = myAssembler.getScalar(energy) // (2)
```

1. This line represents the construction of the SparseFlatAssembler as explained above.
2. To learn what alternatives for `energy` are available and how this works, read the [FE requirements page](feRequirements.md).


## VectorFlatAssembler
It offers the functions of [ScalarAssembler](#scalarassembler) plus additionally
```cpp
Eigen::VectorXd& getVector(const RequirementType& fErequirements)
Eigen::VectorXd& getReducedVector(const RequirementType& fErequirements)
```
As the name suggests, you can either get the full vector or the reduced vector where boundary conditions are considered.
They work the same way as the scalar assembling functions of [ScalarAssembler](#scalarassembler).
The available requirements are explained on the [FE requirements page](feRequirements.md).


## SparseFlatAssembler
It offers the functions of [VectorFlatAssembler](#vectorflatassembler) plus additionally
```cpp
Eigen::SparseMatrix<double> &getMatrix(const RequirementType &fErequirements)
Eigen::SparseMatrix<double> &getReducedMatrix(const RequirementType &fErequirements)
```
A sparse matrix is returned.
They work the same way as the vector assembling functions of [VectorFlatAssembler](#vectorflatassembler).
The available requirements are explained on the [FE requirements page](feRequirements.md).



## DenseFlatAssembler
The only difference between the [SparseFlatAssembler](#sparseflatassembler) and the DenseFlatAssembler is that the
DenseFlatAssembler returns a dense matrix.
```cpp
Eigen::MatrixXd &getMatrix(const RequirementType &fErequirements)
Eigen::MatrixXd &getReducedMatrix(const RequirementType &fErequirements)
```