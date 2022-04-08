# Solvers

In Ikarus there are essentially two types of solvers.

## Linear solver
The first are called `LinearSolver`. These are solver which solve for the vector \( \boldsymbol{x} \) in

$$
\boldsymbol{A}  \boldsymbol{x} =  \boldsymbol{b}
$$
where \(\boldsymbol{A} \) is some matrix and \(\boldsymbol{b}\) is some vector.
These solvers can be direct or iterative. Furthermore, they depend on the underlying structure of the matrix \(\boldsymbol{A} \).
I.e. if it is stored in `dense` or in a `sparse` format.

Currently, we only support the linear solvers provided by the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library.

Linear solvers can be constructed by calling the constructor

```cpp
ILinearSolver(const SolverTypeTag& solverTypeTag)
```

There exits an enum type `SolverTypeTag` with the following values

{{ inputcpp('src/include/ikarus/solver/linearSolver/linearSolver.hh',0,20,40) }}

The prefixes `s_` and  `d_` indicate wether the linear solver can be used for dense or sparse matrices.
Furthermore, there is a second prefix for sparse solvers[^1] `d` and `i` for direct solvers and for iterative solvers.
Thus, using `si_ConjugateGradient` means that this solver is for sparse matrices and is an iterative solver.

[^1]: Dense solver are currently all direct solvers. Therefore, we do not distinguish them.

The naming of the solvers is the same as in Eigen. 
For more details on the solvers we refer to [Eigen's documentation for dense decompositions](https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html)
and to [Eigen's documentation for sparse decompositions](https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html).

### Interface 
Similar to Eigen's interface the following function are provided
```cpp
void analyzePattern(const MatrixType& A);  // (1)
void factorize(const MatrixType& A); // (2)
ILinearSolver& compute(const MatrixType& A); // (3)
void solve(Eigen::VectorX<ScalarType>&x, const Eigen::VectorX<ScalarType>& b); // (4)
```
1. If the matrix is sparse Eigen can collect information on the sparsity pattern of the matrix for faster a faster `solve` step. This pattrern does not change if you change the values of the non-zero entries.
2. This method applies some decomposition for direct solvers e.g. LU decomposition. For iterative solvers the method is a noOp.
3. Compute simply calls 'analyzePattern' and 'factorize'.
4. Solves the problem and stores the result in `x`.

!!! note 

    If your algorithm in mind does rely on special features of some linear solver then you have to directly use this solver.
    E.g. if you need the `.determinant()` method of `Eigen::SimplicialLDLT` you need to directly use it since `ILinearSolver`does not support this method.


## Non-Linear solver

Non-linear solvers are usually used to solve some optimization problem, e.g. root-finding or minimization problems
\begin{align}
  \boldsymbol{R}(\boldsymbol{x}) \stackrel{!}{=} \boldsymbol{0} \quad \text{or} \quad
  \min_{\boldsymbol{x} \in \mathcal{M}} f(\boldsymbol{x} )
\end{align}

### Interface
```cpp
void setup(const NewtonRaphsonSettings& p_settings); // (1)
SolverInformation solve(const SolutionType& dx_predictor = NoPredictor{}); // (2)
auto& nonLinearOperator(); // (3)
```

1. With this function several properties of the nonlinear solver can be set. E.g. residual tolerance or maximum number of iterations.
2. Solves the non-linear problem. One can pass an initial guess to the function. Otherwise the zero vector is assumed. 
   It returns`SolverInformation` which contains information on the sucess of the solution step and other information as the needed iterations. 
3. Simply returns the underlying 

!!! note
    To easy the construction process the Nonlinear solver can provide a method `make[...]` which allows shorter syntax.
    since no `std::shared_ptr` has to be constructed and specifying all template arguments. 
    The construction of the nonlinear solvers can be very differnt therefore we do not impose and interface for the constructors.
### Implementations
| Name                      | Purpose         | Constraints on nonlinear operator                       | Header |Properties |
|:--------------------------|:----------------|:----------------|--|--|
| Newton-Raphson            | Root finding    | Value and gradient          | `newtonRaphson.hh`| Locally quadratic convergence |
| Trust-Region              | Minimization    | Value, gradient and Hessian    | `trustRegion.hh`| Globally convergent and locally quadratic convergence |

To see the Newton-Raphson we refer to the tests inside `nonLinearOperatorTest.cpp` and for trust region `trustRegionTest.cpp`.

\bibliography