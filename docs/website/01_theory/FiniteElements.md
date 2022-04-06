# Finite elements

Several disciplines associate to finite elements different meanings.
In Ikarus finite elements have two different tasks.
The first one is to provide the evaluation of scalars, vectors and matrices. 
These are associated to an algebraic representation of discrete energies, weak forms or bilinear forms
These algebraic objects are usually constructed using some combination of [local function](LocalFunctions.md) and 
parameters steeming from the underlying physical problem, e.g. load factor, Young's modulus or viscocity.

The second task of finite elements is to evaluate derived results in the element parameter space. E.g. stresses or geometric quantities.

## Interface
Local functions provide the following interface
```cpp
void evaluateScalar(const DomainType& local);
void evaluateVector(const unsigned int& integrationPointIndex);
void evaluateMatrix(const DomainType& local,...);
void globalIndices(const unsigned int& integrationPointIndex,...);
auto calculateAt();
```

\bibliography