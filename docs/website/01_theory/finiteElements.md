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

## Linear elasticity
* To be added

## Enhanced Assumed Strain Elements
The Enhanced Assumed Strain (EAS) elements are a class of finite elements which helps to avoid the locking phenomena.
They are obtained by re-parametrizing the Hu-Washizu principle and enforcing an orthogonality condition. 
This results in an extension of the standard pure displacement formulation with an enhanced strain field ($\tilde\epsilon$) 
as an additional independent variable. With an appropriate choice of ansatz space for $\tilde\epsilon$, the locking 
characteristics of the pure displacement formulations can be eliminated. For further theoretical aspects, the readers are referred to [@simo_class_1990] 
and [@andelfinger_eas-elements_1993]. The EAS formulation is currently implemented for the linear elastic case, but 
it could be extended to the non-linear regime. The currently implemented EAS elements are the following:

* Q1E4
* Q1E5
* Q1E7
* H1E9
* H1E21

The notation used here is described as follows. The first alphabet stands for a Quadrilateral (Q) or a Hexahedral (H) element.
The second index denotes the order of the element. E stands for the EAS element and the number following that denotes the 
number of EAS parameters used to enhance the strain field. The only difference amongst various EAS formulations arises 
from the matrix $\mathbf{M}$ which is used to approximate the enhanced strain field. An example for the calculation of the 
matrix $\mathbf{M}$ for a Q1E4 element is shown below:
```cpp
template <typename Geometry>
struct EASQ1E4 {
  static constexpr int strainSize         = 3;
  static constexpr int enhancedStrainSize = 4;

  EASQ1E4() = default;
  explicit EASQ1E4(const Geometry& geometry)
      : geometry{std::make_unique<Geometry>(geometry)}, T0InverseTransformed{calcTransformationMatrix2D(geometry)} {}

  auto calcM(const Dune::FieldVector<double, 2>& quadPos) const {
    Eigen::Matrix<double, strainSize, enhancedStrainSize> M;
    M.setZero(strainSize, enhancedStrainSize);
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * xi - 1.0;
    M(2, 3)           = 2 * eta - 1.0;
    const double detJ = geometry->integrationElement(quadPos);
    M                 = T0InverseTransformed / detJ * M;
    return M;
  }

  std::unique_ptr<Geometry> geometry;
  Eigen::Matrix3d T0InverseTransformed;
};
```
It is to note that the ansatz spaces for the matrix $\mathbf{M}$ are to be modified such that it fulfills the orthogonality 
condition in the $\left[0,1\right]$ element domain used in DUNE, in contrast to the $\left[-1,1\right]$ usually found in 
literature.

In order to add a new EAS element, 

\bibliography