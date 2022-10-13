# Finite elements

Several disciplines associate to finite elements different meanings.
In Ikarus finite elements have two different tasks.
The first one is to provide the evaluation of scalars, vectors and matrices. 
These are associated to an algebraic representation of discrete energies, weak forms or bilinear forms
These algebraic objects are usually constructed using some combination of [local function](localFunctions.md) and 
parameters steeming from the underlying physical problem, e.g. load factor, Young's modulus or viscocity.

The second task of finite elements is to evaluate derived results in the element parameter space. E.g. stresses or geometric quantities.
This boils down to the following interface.
## Interface
Local functions provide the following interface
```cpp
ScalarType evaluateScalar(const FErequirements& req);
void evaluateVector(const FErequirements& req, VectorType& b);
void evaluateMatrix(const FErequirements& req, MatrixType& A);
void calculateLocalSystem(const FErequirements& req, MatrixType& A, VectorType& b);
void calculateAt(const Resultrequirements& req, const Eigen::Vector<double, Traits::mydim>& local,
                     ResultTypeMap<ScalarType>& result);
void globalIndices(std::vector<GlobalIndex>& globalIndices);
```

To discuss these methods first finite element requirements and result requirements should be learned, see [fe requirements](feRequirements.md).
The first four methods receive an object of type `FErequirements`. This object is responsible for passing different information needed for the local evaluation of the local linear algebra objects.
The first method `evaluateScalar` simply returns by values since usually this is cheap to return a `double`.
The other methods `evaluateVector`, `evaluateMatrix` and `calculateLocalSystem` receive one or two  output argument where the result should be written.
This interface is needed to circumvent the dynamic memory allocation, if these methods would return by value.

The method `calculateAt` is responable to evaluate several results and it receives a `ResultRequirements` object which contains information which results should be evaluated.
These results are stored inside the output argument `result` which is of type `ResultTypeMap`.
Additionally there is the argument 'local' which stores the coordinates where inside the element coordinates the result should be evaluated.


Inside a typical `calculateAt` method the usage is 

```cpp
typename ResultTypeMap<double>::ResultArray res;
if(req.isResultRequested( ResultType::gradientNormOfMagnetization)) {
  res.resize(1,1);
  res(0,0)=...;
  result.insertOrAssignResult(ResultType::gradientNormOfMagnetization,res);
}
if(req.isResultRequested( ResultType::BField)) {
  res.setZero(3,1);
  res=...;
  result.insertOrAssignResult(ResultType::BField,res);
}
if(req.isResultRequested( ResultType::cauchyStress)) {
  res.setZero(3,3);
  res = ...;
  result.insertOrAssignResult(ResultType::cauchyStress,res);
}
```
!!! note "`ResultTypeMap<double>::ResultArray`"
    `#!cpp ResultTypeMap<double>::ResultArray` is an object of type `#!cpp Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,3,3>`.
    Thus, the maximum result size is limited to a 3x3 matrix. This is used to circumvent dynamic memory allocations.


The last method is `globalIndices`. It is used to message the global indices of this finite element steming in the output parameter `globalIndices`.
This information should stem from a basis object. See existing implementations for details.

\bibliography