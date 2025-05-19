# Solvers

In Ikarus, there are essentially two types of solvers: linear and non-linear.

## Linear solver

It solves for the vector \( \boldsymbol{x} \) in

$$
\boldsymbol{A}  \boldsymbol{x} =  \boldsymbol{b}
$$
where \(\boldsymbol{A} \) and \(\boldsymbol{b}\) are any matrix or vector, respectively.
These solvers can be direct or iterative. Furthermore, they depend on the underlying structure of the
matrix \(\boldsymbol{A} \), i.e., whether it is stored in a `dense` or `sparse` format.

Currently, only the linear solvers provided by the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
library are supported.

Linear solvers can be constructed by calling the constructor:

```cpp
LinearSolver(const SolverTypeTag& solverTypeTag)
```

There exists an enum type `SolverTypeTag` with the following values:

{{ inputcpp('ikarus/solver/linearsolver/linearsolver.hh',False,19,39) }}

The prefixes `s_` and  `d_` indicate whether the linear solver can be used for `sparse` or `dense` matrices, respectively.
Furthermore, there is also a second prefix for sparse solvers: `d` and `i` for direct solvers and for iterative solvers, respectively.
Thus, using `si_ConjugateGradient` means that this solver is for sparse matrices and is an iterative solver.

!!! note
    All dense solvers are currently direct solvers. Therefore, we do not distinguish them.

The naming of the solvers is the same as in the Eigen library.
For more details on the solvers, we refer to [Eigen's documentation for dense decompositions](https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html)
and to [Eigen's documentation for sparse decompositions](https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html).

### Interface

Similar to Eigen's interface, the following functions are provided:

```cpp
void analyzePattern(const MatrixType& A);  // (1)!
void factorize(const MatrixType& A); // (2)!
LinearSolver& compute(const MatrixType& A); // (3)!
void solve(Eigen::VectorX<ScalarType>&x, const Eigen::VectorX<ScalarType>& b); // (4)!
```

1. If the matrix is sparse, Eigen can collect information on the sparsity pattern of the matrix for a faster `solve` step. This pattern
   does not change if the non-zero entries are modified.
2. This method applies a decomposition technique for direct solvers, e.g., LU decomposition. For iterative solvers, the method is a
   non-linear operator.
3. It calls both the functions `analyzePattern` and `factorize`.
4. Solves the problem and stores the result in `x`.

!!! tip
    If your algorithm relies on special features or attributes of a linear solver, then the solver is to be directly used.
    For example, if the `.determinant()` method of `Eigen::SimplicialLDLT` is required, it must be called directly because
    `LinearSolver`does not support it.

## Nonlinear solver

Non-linear solvers are usually used to solve a kind of optimization problem, e.g., root-finding or minimization problems like:
\begin{align}
  \boldsymbol{R}(\boldsymbol{x}) \stackrel{!}{=} \boldsymbol{0} \quad \text{or} \quad
  \min_{\boldsymbol{x} \in \mathcal{M}} f(\boldsymbol{x} )
\end{align}

It has the following interface:

```cpp
void setup(const NewtonRaphsonSettings& p_settings); // (1)!
SolverInformation solve(const SolutionType& dx_predictor = NoPredictor{}); // (2)!
auto& differentiableFunction(); // (3)!
```

1. With this function, several properties of the nonlinear solver can be set. E.g., residual tolerance or maximum number of iterations
2. Solves the non-linear problem. An initial guess to the function can be passed, otherwise, a zero vector is assumed.
   It returns the `SolverInformation` which contains information like the success of the solution step and more.
3. Just returns the underlying `differentiableFunction`, see [link](differentiablefunction.md).

!!! note
    To ease the construction process, the non-linear solver can provide a method `make[...]` that allows shorter syntax
    since no `std::shared_ptr` has to be constructed specifying the template arguments.
    The construction of the nonlinear solvers can be very different. Therefore, we do not impose an interface for the constructors.

### Implementations

| Name                      | Purpose         | Constraints on the function                       | Header                                         |Properties |
|:--------------------------|:----------------|:----------------|------------------------------------------------|--|
| Newton-Raphson            | Root finding    | Value and gradient          | `newtonRaphson.hh`                             | Locally quadratic convergence |
| Newton-Raphson with scalar subsidiary function            | Root finding with a scalar function as additional constraint    | Value, gradient and a scalar function          | `newtonraphsonwithscalarsubsidiaryfunction.hh` | Locally quadratic convergence |
| Trust-Region              | Minimization    | Value, gradient and Hessian    | `trustRegion.hh`                               | Globally convergent and locally quadratic convergence |

To see the Newton-Raphson implementation, we refer to the tests inside `differentiableFunctionTest.cpp`, and for the
trust-region method, `trustRegionTest.cpp`.

\bibliography
