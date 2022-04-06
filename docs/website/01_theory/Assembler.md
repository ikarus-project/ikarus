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
- `dirichFlags` is of type `std::vector<bool>`. `dirichFlags[i] = 1` means that degree of freedom i is fixed. 
    The corresponding row / column / entry will be eliminated when you ask for reduced matrix / vector. 

## SparseFlatAssembler
The SparseFlatAssembler has to public member functions:
```cpp
Eigen::SparseMatrix<double>& getMatrix(const RequirementType& fErequirements)
```
It assembles the reqested matrix and returns a sparse matrix. A call to this function could look as follows:
```cpp
SparseFlatAssembler myAssembler(...) // (1)
// other code
const auto& K = myAssembler.getMatrix(stiffness) // (2)
```

1. This line represents the construction of the SparseFlatAssembler as explained above.
2. To learn what alternatives for `stiffness` are available and how this works, read [the available requirements section](#available-requirements) below.

```cpp
Eigen::SparseMatrix<double>& getReducedMatrix(const RequirementType& fErequirements)
```
returns the reduced matrix, i.e. the rows and columns which are flagged with `1` in `dirichFlags` are eliminated. 

## VectorFlatAssembler
It offers the functions
```cpp
Eigen::VectorXd& getVector(const RequirementType& fErequirements)
Eigen::VectorXd& getReducedVector(const RequirementType& fErequirements)
```

They work the same way as the matrix assembling functions of [SparseFlatAssembler](#sparseflatassembler).
The available requirements are explained in [the available requirements section](#available-requirements) below.

## ScalarAssembler
It offers only one function
```cpp
double& getScalar(const RequirementType& fErequirements)
```
The available requirements are explained in [the available requirements section](#available-requirements) below.
`dirichletFlags` is not used in this assembler.

## DenseFlatSimpleAssembler

## DenseFlatAssembler

## Available requirements