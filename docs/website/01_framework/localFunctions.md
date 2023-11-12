# Local functions

Local functions are functions that are bound to single grid elements.
Therefore, they are constructed from some local basis, a coefficient vector, and the geometry of the grid element.
Since the implementation is quite involved, `localfefunctions` do not reside at Ikarus but in the separate
module [dune-localfefunctions](https://github.com/ikarus-project/dune-localfefunctions).

Usually, local functions are need to be evaluated in the local coordinate system \( \mathbb{\xi} \in T_{\text{ref}} \subset\mathbb{R}^n \) :

$$
f: \boldsymbol{\xi}^n \rightarrow \mathbb{R}^m
$$

where $T_{\text{ref}}$ is the reference element, e.g., for a hypercube, $T_{\text{ref}}= [0,1]^d$.

## Interface

Local functions provide the following interface:

```cpp
LocalFunction(const Dune::CachedLocalBasis<DuneBasis>& p_basis, const CoeffContainer& coeffs_,
              const std::shared_ptr<const Geometry>& geo,
              Dune::template index_constant<ID> = Dune::template index_constant<std::size_t(0)>{}); // (1)!

FunctionReturnType evaluate(const DomainType& local);
FunctionReturnType evaluate(const unsigned int& integrationPointIndex);

auto evaluateDerivative(const DomainType& local,...);
auto evaluateDerivative(const unsigned int& integrationPointIndex,...);
auto viewOverIntegrationPoints(); // (2)!

template<std::size_t ID=0>
constexpr int order(Dune::index_constant<ID> ); // (3)!

template<std::size_t ID=0>
auto basis(Dune::index_constant<ID> ); // (4)!

template<std::size_t ID=0>
auto coefficientsRef(Dune::index_constant<ID>); // (5)!

template <typename IntegrationRule, typename... Ints>
void bind(IntegrationRule&& p_rule, Derivatives<Ints...>&& ints); // (6)!

auto clone (); // (7)!

template<typename ScalarType, std::size_t ID=0>
auto rebindClone (ScalarType, Dune::index_constant<ID>); // (8)!
```

1. The constructor takes a `Dune::CachedLocalBasis`, a vector of coefficients and a shared pointer to the geometry of the grid elements.
   Additionally, the local function can be tagged with a compile-time constant, e.g., `Dune::template index_constant<0>` or
   `Dune::Indices::_0`.
2. This returns a vector of structs representing the integration point and its index. Therefore, the syntax is `#!cpp for (const auto&
   [gpIndex, gp] : localFunction.viewOverIntegrationPoints()) {...}`
3. Returns the order of the local function with respect to the coefficients. An id tag can be passed, which returns the order with
   respect to a tagged function. For details, see [Tagging leaf local functions](#tagging-leaf-local-functions).
4. Returns the basis of the local function. An id tag can be passed, which returns the basis of a specific tagged function. For details
   see [Tagging leaf local functions](#tagging-leaf-local-functions).
5. Returns a reference to the coefficient of the underlying leaf local finite elements. An id tag can be passed, which returns the basis
   of a specific tagged function. It can return a const or non-const reference. The non-const version is deactivated, if there is more
   than one leaf node with the passed id tag.  For details, see [Tagging leaf local functions](#tagging-leaf-local-functions).
6. This function is passed to the given `localBasis`. See [link](localBasis.md)
7. Clones the local function and stores a copy of all leave nodes.
8. Clones the local function and rebinds the scalar type of the coefficients with id tag `ID`. This comes in handy if, for example, one
   wants to replace doubles with an autodiff type.

The "..." in the `evaluateDerivative` function call refers to several possible variadic templates.
The implementation looks like the following:

=== "Usage with integration point index"

    ```cpp
    using namespace Dune::DerivativeDirections;
    localFunction.bind(rule, bindDerivatives(0,1));
    for(auto& [gpIndex, gp] : localFunction.viewOverIntegrationPoints()){
      localFunction.evaluateDerivative(gpIndex, wrt(spatialAll)); // (1)!
      localFunction.evaluateDerivative(gpIndex, wrt(spatialAll), on(gridElement)); // (2)!
    }
    ```

    1. Compute the spatial Jacobian of localFunction
    2. Compute the spatial Jacobian of localFunction and transform it to the grid element

=== "using integration point coordinates"

    ```cpp
    using namespace Dune::DerivativeDirections;
    for(auto& gp : rule){
      localFunction.evaluateDerivative(gp.position(), wrt(spatialAll)); // (1)!
      localFunction.evaluateDerivative(gp.position(), wrt(spatialAll), on(gridElement)); // (2)!
    }
    ```

    1. Compute the spatial Jacobian of localFunction
    2. Compute the spatial Jacobian of localFunction and transform it to the grid element

where the first call implements

$$
\operatorname{grad}_\boldsymbol{\xi} f : \boldsymbol{\xi} \rightarrow \mathbb{R}^{m \times d}.
$$

The second one takes into account the fact that the local function is defined in some physical space $\boldsymbol{X}$ with the
coordinate $\boldsymbol{x}$.
Therefore, it transforms the Jacobian from the reference element $\operatorname{grad}_{\boldsymbol{\xi}}$ to the Jacobian on the grid
element $\operatorname{grad}_\boldsymbol{x}$.
This behavior is activated, if `on(gridElement)` is passed; otherwise, if the derivatives on the reference element is needed, pass `on
(referenceElement)`.
Here, `gridElement` and `referenceElement` are global constants in the namespace `Dune::DerivativeDirections`.

Thus, if `on(gridElement)` is passed, the local function usually implements

$$
\operatorname{grad}_\boldsymbol{x} = \operatorname{grad}_{\boldsymbol{\xi}} \boldsymbol{J}^{-1}
$$

where $J$ is the Jacobian of the mapping from the reference element $T_{\text{ref}}$ to the element living in physical space $T$.
For details, see page 22 of the Dune book[@sander2020dune].

Instead of passing `spatialAll` to `wrt(..)`, there are other helper functions such as:

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(spatial(0))); // (1)!
localFunction.evaluateDerivative(gpIndex, wrt(spatial(1))); // (2)!
```

1. Compute the first column of the spatial Jacobian of `localFunction`
2. Compute the second column of the spatial Jacobian of `localFunction`

which can also be combined with `on(...)`.

## Derivatives w.r.t. coefficients

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(coeff(j)));
```

which evaluates the first derivative for a vector space valued function,e.g., for $f(\boldsymbol{\xi}) = \sum_{I=1}^n N^I(\boldsymbol{\xi}) \boldsymbol{x}_I$,
we arrive at a matrix $\boldsymbol{A}$ such that

$$
[\boldsymbol{A}]_{ij}  = A_{ij} =  \frac{\partial f_i(\boldsymbol{\xi})}{\partial \boldsymbol{x}_j}
$$

and the second derivative

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k)), along(q));
```

$$
[\boldsymbol{B}]_{jk} =  B_{jk} = q_i A_{ijk} =  \frac{\partial^2 (q_i  f_i(\boldsymbol{\xi}))}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}
$$

where $\boldsymbol{q}$ is an arbitrary vector of the same size as $f$, i.e., it is the direction of the derivative in this case.
$\boldsymbol{A}$ and $\boldsymbol{B}$ are simply then the returned matrices and do not have any special meaning.
If a vector is not passed while evaluating the second derivative, the result would be a third-order tensor for a vector-valued function $f$.
As a result, a direction derivative in the direction given by `along(q)` is computed to return a matrix $\boldsymbol{B}$ in this case.
This helps for both readability and performance. See the [example](#example-dirichlet-energy) later for more details.

## Derivatives w.r.t. coefficients and spatial derivatives

Spatial derivatives and derivatives w.r.t. the coefficients can be combined. Therefore, it is legal to call

```cpp
auto B = localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k),spatialAll), along(Q));
auto b1 = localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k),spatial(0)), along(q));
auto b2 = localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k),spatial(1)), along(q));
```

!!! warning

    The order of spatial and coefficient derivatives does not matter. The returned value is always rearranged so that the first derivative is the spatial one.

The first line is then equivalent to

$$
[\boldsymbol{B}]_{jk} =  B_{jk} = Q_{il} A_{iljk} =  \frac{\partial^2
([\operatorname{grad}_\boldsymbol{\xi} f(\boldsymbol{\xi})]_{il} Q_{il} )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
$$

For the second and third line, we have

\begin{align}
\boldsymbol{b}_{0,jk} = \frac{\partial^2
([\operatorname{grad}_{\xi^0} f(\xi)]_{i} q_i )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}, \\
\boldsymbol{b}_{1,jk} = \frac{\partial^2 ([\operatorname{grad}_{\xi^1} f(\xi)]_{i} q_i )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
\end{align}

These objects are also returned when the second and third lines above are used.

All of these function calls, once again, can be combined with `on(gridElement)` as shown below:

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k),spatialAll), along(Q), on(gridElement));
```

which computes

$$
\frac{\partial^2 ([\operatorname{grad}_\boldsymbol{x} f(\boldsymbol{\xi})]_{il} Q_{il} )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
$$

!!! warning "Warning"

    Currently, only first-order spatial derivatives and second-order derivatives w.r.t. the coefficients are supported.

## Example: Dirichlet energy

This example shows how the energy, gradient, and Hessian of a [dirichlet energy](https://en.wikipedia.org/wiki/Dirichlet_energy) can be calculated.
$$
E(\boldsymbol{u}) = \frac{1}{2} \int_\Omega ||\operatorname{grad}_\boldsymbol{x} \boldsymbol{u}(\boldsymbol{x})|| ^2 \textrm{d} \boldsymbol{x}
$$

If the energy is to be minimized w.r.t. the coefficients of the nodes, the energy, gradient, and Hessian w.r.t. the coefficients are to be calculated.
Of course, this depends on the optimization algorithms, but for now, the general case where all three are needed is considered.

```cpp
auto dirichletEnergy() {
  double energy = 0;
  // bind localBasis to some integration rule
  // create uNodalCoeffs
  Ikarus::StandardLocalFunction uFunction(localBasis, uNodalCoeffs, sharedGeometry);
  for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
    //.. calculate the inverse Jacobian of the geometry
    const auto gradu = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll), on(gridElement));
    energy+= 0.5 * (gradu.transpose() * gradu).trace() * gp.weight() * sharedGeometry->integrationElement(gp.position());
  }
}
```

```cpp
auto gradientDirichletEnergy(Eigen::VectorXd& g) {
  //.bind localBasis to some integration rule
  // create uNodalCoeffs
  constexpr int size =  // spatial size of u
      Ikarus::StandardLocalFunction uFunction(localBasis, uNodalCoeffs, sharedGeometry);
  for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
    //.. calculate the inverse Jacobian of the geometry
    const auto gradu = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll), on(gridElement));
    for (auto i : fe.size()) { //loop over coeffs, i.e.nodes of the finite element
      const auto graduDCoeffs
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)), on(gridElement));
      Eigen::Vector<double, size> tmp;
      tmp.setZero();
      for (int k = 0; k < gridDimension; ++k)
        tmp += graduDCoeffs[k] * gradu.col(k);  // (1)!
      g.segment<size>(i * size) += tmp * gp.weight() * sharedGeometry->integrationElement(gp.position());
    }
  }
}
```

1. `graduDCoeffs` contains in `graduDCoeffs[0]` the derivatives w.r.t. the coefficient of the first column, and at `graduDCoeffs[1]` the
   derivatives w.r.t. the second column of `gradu`.

```cpp
auto hessianDirichletEnergy(Matrix& h) {
  //... bind localBasis to some integration rule
  // and create uNodalCoeffs
  constexpr int size =  // spatial size of u
      Ikarus::StandardLocalFunction uFunction(localBasis, uNodalCoeffs, sharedGeometry);
  for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
    //.. calculate the inverse Jacobian of the geometry
    for (auto i : loop over coeffs, i.e.nodes of the finite element) {
      const auto graduDCoeffsI
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)), on(gridElement));
      for (auto j : fe.size()) { //loop over coeffs, i.e.nodes of the finite element
        const auto graduDCoeffsJ
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeffs), on(gridElement), coeffIndices(j));
        Eigen::Matrix<double, size, size> tmp;
        tmp.setZero();
        for (int k = 0; k < gridDimension; ++k)
          tmp += graduDCoeffsI[k] * graduDCoeffsJ[k];
        h.block<size, size>(i * size, j * size) += tmp * gp.weight() * sharedGeometry->integrationElement(gp.position());
      }
    }
  }
}
```

## Implementations

In the following, the local functions that are currently available are summarized.
The ansatz functions are denoted as $N^i(\boldsymbol{\xi})$ in the table below.

| Name                                               | Interpolation formula                                                                     | Note                                                                                                                                                                                                                                                         | Header                            |
|:---------------------------------------------------|:------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------|
| Standard                                           | $$ \boldsymbol{x} = \sum_{i=1}^n N^i(\boldsymbol{\xi}) \boldsymbol{x}_i  $$               | -                                                                                                                                                                                                                                                            | `standardLocalFunction.hh`        |
| Projection-Based[@grohs_ProjectionBasedFinite2019] | $$ \boldsymbol{x} = P\left(\sum_{i=1}^n N^i(\boldsymbol{\xi}) \boldsymbol{x}_i \right) $$ | This is one version of geometric finite elements. These are finite elements suited <br /> for interpolation on manifolds. Here $P: \mathbb{R}^m \rightarrow \mathcal{M}$ is an operator that projects <br /> the usual linear interpolation onto a manifold. | `projectionBasedLocalFunction.hh` |

## How to implement your own local functions

To implement your own local function, the file
[`ikarus/localFunctions/impl/localFunctionTemplate.hh`](https://github.com/IkarusRepo/Ikarus/src/include/ikarus/localFunctions/impl/localFunctionTemplate.hh)
is made available.

The file can be copied, and then you can rename the class to your preferred name and implement the functions mentioned below.
If a particular function is not required, it has to be deleted explicitly.
Then, if someone calls that function, it returns a zero matrix or vector of appropriate size.

These functions are all templated with `DomainTypeOrIntegrationPointIndex`, which is an integration point index or position.
Additionally, `On<TransformArgs>` specifies whether the function should be evaluated on the reference element or the grid element (see
above).

```cpp
FunctionReturnType evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                        const On<TransformArgs>&) const; // (1)!

Jacobian evaluateDerivativeWRTSpaceAllImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                           const On<TransformArgs>& transArgs) const; // (2)!

JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                     int spaceIndex, const On<TransformArgs>& transArgs) const; // (3)!


CoeffDerivMatrix evaluateDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                 int coeffsIndex, const On<TransformArgs>& transArgs) const; // (4)!

CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                       const std::array<size_t, 2>& coeffsIndex,
                                                       const Along<AlongArgs...>& alongArgs,
                                                       const On<TransformArgs>& transArgs) const; // (5)!

std::array<CoeffDerivEukRieMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
    const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
    const On<TransformArgs>& transArgs) const; // (6)!


CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
    const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
    const On<TransformArgs>& transArgs) const;  // (7)!


auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
    const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
    const Along<AlongArgs...>& alongArgs, const On<TransformArgs>& transArgs) const; // (8)!

CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
    const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
    const int spatialIndex, const Along<AlongArgs...>& alongArgs, const On<TransformArgs>& transArgs) const; // (9)!
```

1. This is called by `localFunction.evaluate(...)`.
2. This is called by `localFunction.evaluateDerivative(..., wrt(spatialAll))`.
3. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i)))`.
4. This is called by `localFunction.evaluateDerivative(..., wrt(coeff(j)))`.
5. This is called by `localFunction.evaluateDerivative(..., wrt(coeff(j,k)))`.
6. This is called by `localFunction.evaluateDerivative(..., wrt(spatialAll,coeff(j)))`.
7. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i),coeff(j)))`.
8. This is called by `localFunction.evaluateDerivative(..., wrt(spatialAll,coeff(j,k)), along(A))`. `A` can be accessed via `std::get<0>(alongArgs.args)`.
9. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i),coeff(j,k)), along(v))`. `v` can be accessed via `std::get<0>(alongArgs.args)`.

## Expressions

[Expression templates](https://de.wikipedia.org/wiki/Expression_Templates) are usually used in linear algebra libraries, e.g.,
[Eigen](https://eigen.tuxfamily.org) or [Blaze](https://bitbucket.org/blaze-lib/blaze/src/master/).
The syntax is similar to the one provided by [UML](https://fenics.readthedocs.io/projects/ufl/en/latest/manual/form_language.html) but
only acts on local functions.
Expression templates are used here to combine existing local functions in order to obtain new nested ones.

For example, consider the following code:

```cpp
...
auto f = Ikarus::StandardLocalFunction(localBasis, coeffVectors0, sharedGeometry);
auto g = Ikarus::StandardLocalFunction(localBasis, coeffVectors1, sharedGeometry);
```

Two local functions that satisfy the interface described above are created.
Now it is possible to combine these functions and get an object that also satisfies the concept above.
Thus, the following is possible:

```cpp
...
auto f = Ikarus::StandardLocalFunction(localBasis, coeffVectors0, sharedGeometry);
auto g = Ikarus::StandardLocalFunction(localBasis, coeffVectors1, sharedGeometry);
auto k = f + g;
k.evaluateDerivative(ipIndex, wrt(coeff(i), spatial(d)));
```

Currently, binary and unary expressions are supported. The following expressions are defined:

| Name                   | Mathematical formula                                                                                                                                 | Code                              | Note                                                                                                    |
|:-----------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------|:--------------------------------------------------------------------------------------------------------|
| Sum                    | $$ \boldsymbol{f} + \boldsymbol{g}  $$                                                                                                               | `#!cpp f+g`                       | $\boldsymbol{f}$  and  $\boldsymbol{g}$ needs to be of the same size.                                   |
| DotProduct             | $$ \boldsymbol{f} \cdot \boldsymbol{g} = f_i g_i $$                                                                                                  | `#!cpp dot(f,g)`                  | $\boldsymbol{f}$  and  $\boldsymbol{g}$ needs to be of the same size.                                   |
| normSquared            | $$ \boldsymbol{f} \cdot \boldsymbol{f} = f_i f_i $$                                                                                                  | `#!cpp normSquared(f)`            |                                                                                                         |
| Negate                 | $$ -\boldsymbol{f}  $$                                                                                                                               | `#!cpp -f`                        |                                                                                                         |
| sqrt                   | $$ \sqrt{f}  $$                                                                                                                                      | `#!cpp sqrt(f)`                   | The function $f$ needs a scalar return type.                                                            |
| log                    | $$ \log{f}  $$                                                                                                                                       | `#!log log(f)`                    | The function $f$ needs a scalar return type. Here, log is the natural logarithm.                        |
| pow                    | $$ f^n  $$                                                                                                                                           | `#!cpp pow<n>(f)`                 | The function $f$ needs a scalar return type. $n$ is an integer given during compile-time.               |
| Scale                  | $$  a f , \quad a \in  \mathbf{R}$$                                                                                                                  | `#!cpp a*f` and `#!cpp f/a`       | `#!cpp a` has to satisfy `#!cpp std::is_arithmetic<..>`                                                 |
| LinearStrains          | $$ \frac{1}{2}\left(\boldsymbol{H}+\boldsymbol{H}^T \right),\quad \boldsymbol{H} = \mathrm{grad}(\boldsymbol{f})  $$                                 | `#!cpp linearStrains(f)`          | If you call the formula on the left with `on(gridElement)`, it will assume the transformed derivatives. |
| GreenLagrangianStrains | $$ \frac{1}{2}\left(\boldsymbol{H}+\boldsymbol{H}^T +\boldsymbol{H}^T \boldsymbol{H}\right),\quad \boldsymbol{H} = \mathrm{grad}(\boldsymbol{f})  $$ | `#!cpp greenLagrangianStrains(f)` | If you call the formula on the left with `on(gridElement)`, it will assume the transformed derivatives. |

These expressions can also be nested. As a result, it is valid to write

```cpp
auto f = Ikarus::StandardLocalFunction(localBasis, coeffVectors0);
auto g = Ikarus::StandardLocalFunction(localBasis, coeffVectors1);
auto k = -sqrt(dot(2*f+f,5*g));
k.evaluateDerivative(ipIndex, wrt(coeff(i), spatial(d)));
```

To use these expressions, there are additional exported static types for all expressions.

```cpp
constexpr bool isLeaf; // (1)!
constexpr bool children; // (2)!
```

1. This is true if the underlying expression is one of the above local functions that really contain the coefficients; see [Implementations](#implementations).
2. Returns the number of children (2 for binary expressions and 1 for unary expressions).

!!! note
    To use these expression, simply include the header `#!cpp #include <dune/localfefunctions/expressions.hh>`.

## Tagging leaf local functions

In the context of mixed finite elements, there are usually several local functions that contribute to the energy. These stem from
different local bases.
For example, consider the Q1P0 element, where displacements are interpolated by using the four bilinear ansatz functions and the
element-wise constant pressure field.

Then, to obtain gradients and Hessians, we need to differentiate w.r.t. different coefficients.
This can be done by tagging the local function during construction.

```cpp
using namespace Dune::Indices;
auto f = Ikarus::StandardLocalFunction(localBasis0, coeffVectors0, sharedGeometry, _0);
auto g = Ikarus::StandardLocalFunction(localBasis1, coeffVectors1, sharedGeometry, _1);
auto k = dot(f,g);
k.evaluateDerivative(ipIndex, wrt(coeff(_0,i,_1,j))); // Second derivative w.r.t. the nodal coefficients
```

To explain the last line above, let us consider that the function f is constructed as $f= \sum_{I=0}^n N^L f_L$ and similarly
$g= \sum_{I=0}^m M^K g_K$, where $N$ and $M$ are ansatz functions and $f_L$ and $g_K$ are nodal coefficients.

Thus, the above call translates to

\begin{align}
\frac{\partial^2 (\boldsymbol{f} \cdot \boldsymbol{g} )}{\partial \boldsymbol{f}_i\partial \boldsymbol{g}_j},
\end{align}

where the correct sizes of the result are derived at compile-time.
The complete Hessian of $\boldsymbol{f} \cdot \boldsymbol{g}$ can be calculated by the following:

```cpp
using namespace Dune::Indices;
auto hessianDirichletEnergy(Matrix& h) {
  //... bind localBasis to some integration rule
  using namespace Dune::Indices;
  auto f = Ikarus::StandardLocalFunction(localBasis0, coeffVectors0, sharedGeometry, _0);
  auto g = Ikarus::StandardLocalFunction(localBasis1, coeffVectors1, sharedGeometry, _1);
  auto k = dot(f,g);
  constexpr int sizef = f.correctionSize; // spatial size of the correction of the coefficients of f
  constexpr int sizeg = g.correctionSize; // spatial size of the correction of the coefficients of g
  constexpr int coeffSizef = coeffVectors0.size();
  constexpr int coeffSizeg = coeffVectors1.size();

  Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<MatrixBlock00,MatrixBlock01>,
                                       Dune::MultiTypeBlockVector<MatrixBlock10,MatrixBlock11> > KBlocked; // (1)!


  for (const auto& [gpIndex, gp] : k.viewOverIntegrationPoints()) {
    for (size_t I = 0; I < coeffSizef; ++I)
      for (size_t J = 0; J < coeffSizef; ++J)
        KBlocked[_0,_0].block<sizef, sizef>(I * sizef, J * sizef)
          += k.evaluateDerivative(gpIndex, wrt(coeff(_0,I,_0,J))) * gp.weight() * sharedGeometry->integrationElement(gp.position());

    for (size_t I = 0; I < coeffSizef; ++I)
      for (size_t J = 0; J < coeffSizeg; ++J)
        KBlocked[_0,_1].block<sizef, sizeg>(I * sizef, J * sizeg)
          += k.evaluateDerivative(gpIndex, wrt(coeff(_0,I,_1,J))) * gp.weight() * sharedGeometry->integrationElement(gp.position());

    for (size_t I = 0; I < coeffSizeg; ++I)
      for (size_t J = 0; J < coeffSizeg; ++J)
        KBlocked[_1,_1].block<sizeg, sizeg>(I * sizeg, J * sizeg)
          += k.evaluateDerivative(gpIndex, wrt(coeff(_1,I,_1,J))) * gp.weight() * sharedGeometry->integrationElement(gp.position());

    for (size_t I = 0; I < coeffSizeg; ++I)
      for (size_t J = 0; J < coeffSizef; ++J)
        KBlocked[_1,_0].block<sizef, sizeg>(I * sizeg, J * sizef)
          += k.evaluateDerivative(gpIndex, wrt(coeff(_1,I,_0,J))) * gp.weight() * sharedGeometry->integrationElement(gp.position());
    }
}
```

1. This block structure is not necessary. Additionally, in this example, all types (MatrixBlock00, MatrixBlock01, MatrixBlock10,
   MatrixBlock11) are considered as `#!cpp Eigen::MatrixXd`.

## Writing your own expressions

It is also possible to write your own expressions. To do so, please take a look at the existing expressions.
The `sqrt` and `normSquared` expressions are the most general unary and binary expressions implemented.

### Implementing the return value

Implementing a return value is the first step in implementing an expression.
This is done by using the following function:

```cpp
template <typename LFArgs>
auto evaluateValueOfExpression(const LFArgs &lfArgs) const;
```

!!! warning
    The interface dictates that the return value needs to be an `Eigen` type. Thus, even if a scalar double is to be returned, it is to
be wrapped in `#!cpp Eigen::Vector<double, 1>`

Additionally, the evaluation of the derivative is to be implemented as shown below:

```cpp
template <int DerivativeOrder, typename LFArgs>
auto evaluateDerivativeOfExpression(const LFArgs &lfArgs) const;
```

### Evaluate the underlying functions

Expressions always act on existing expressions. Therefore, to have the correct return value for the expression, the underlying
quantities are to be evaluated.
`#!cpp this->m()` is used to access unary functions, and `#!cpp this->l()` and `#!cpp this->r()` are used to access binary expressions.

To evaluate unary functions, the following syntax is used:

```cpp
const auto mEvaluated = evaluateFunctionImpl(this->m(), lfArgs);
```

and for binary functions,

```cpp
const auto l_Evaluated = evaluateFunctionImpl(this->l(), lfArgs);
const auto r_Evaluated = evaluateFunctionImpl(this->r(), lfArgs);
```

Because the expression conforms to the syntax of a local function, its derivative can also be evaluated.

In the function `evaluateDerivativeOfExpression`, the template argument `DerivativeOrder` contains the derivative order.
Additionally, the derivative types can also be accessed using the static booleans, as shown below:

```cpp
    static constexpr bool hasTwoCoeff;
    static constexpr bool hasSingleCoeff;
    static constexpr bool hasNoCoeff;
    static constexpr bool hasNoSpatial;
    static constexpr bool hasOneSpatialAll;
    static constexpr bool hasOneSpatialSingle;
    static constexpr bool hasOneSpatial;
```

Using the _dot-product_ as a binary expression example, we have

```cpp
template <int DerivativeOrder, typename LFArgs>
    auto evaluateDerivativeOfExpression(const LFArgs &lfArgs) const {
      const auto u = evaluateFunctionImpl(this->l(), lfArgs);
      const auto v = evaluateFunctionImpl(this->r(), lfArgs);
      if constexpr (DerivativeOrder == 1)  // (1)!
      {
        const auto u_x = evaluateDerivativeImpl(this->l(), lfArgs); // (2)!
        const auto v_x = evaluateDerivativeImpl(this->r(), lfArgs);
        return Ikarus::eval(v.transpose() * u_x + u.transpose() * v_x); // (3)!
      } else if constexpr (DerivativeOrder == 2) {   // (4)!
        const auto &[u_x, u_y] = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs); // (5)!
        const auto &[v_x, v_y] = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs);
        if constexpr (LFArgs::hasNoSpatial and LFArgs::hasTwoCoeff) { // (6)!
          const auto alonguArgs = replaceAlong(lfArgs, along(v)); // (7)!
          const auto alongvArgs = replaceAlong(lfArgs, along(u));

          const auto u_xyAlongv = evaluateDerivativeImpl(this->l(), alongvArgs); // (8)!
          const auto v_xyAlongu = evaluateDerivativeImpl(this->r(), alonguArgs);

          return Ikarus::eval(u_xyAlongv + transpose(u_x) * v_y + transpose(v_x) * u_y + v_xyAlongu);
        } else if constexpr (LFArgs::hasOneSpatial and LFArgs::hasSingleCoeff) { // (9)!
          const auto u_xy = evaluateDerivativeImpl(this->l(), lfArgs); // (10)!
          const auto v_xy = evaluateDerivativeImpl(this->r(), lfArgs);
          if constexpr (LFArgs::hasOneSpatialSingle and LFArgs::hasSingleCoeff) { // (11)!
            return Ikarus::eval(transpose(v) * u_xy + transpose(u_x) * v_y + transpose(v_x) * u_y
                                + transpose(u) * v_xy);
          } else if constexpr (LFArgs::hasOneSpatialAll and LFArgs::hasSingleCoeff) { // (12)!
            std::array<std::remove_cvref_t<decltype(Ikarus::eval(transpose(v) * u_xy[0]))>, gridDim> res; // (13)!
            for (int i = 0; i < gridDim; ++i)
              res[i] = Ikarus::eval(transpose(v) * u_xy[i] + transpose(u_x.col(i)) * v_y + transpose(v_x.col(i)) * u_y
                                    + transpose(u) * v_xy[i]);
            return res;
          }
        }
      } else if constexpr (DerivativeOrder == 3) { // (14)!
        if constexpr (LFArgs::hasOneSpatialSingle) {  // (15)!
          const auto argsForDyz = lfArgs.extractSecondWrtArgOrFirstNonSpatial(); // (16)!

          const auto &[u_x, u_y, u_z] = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs); // (17)!
          const auto &[v_x, v_y, v_z] = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs);
          const auto &[u_xy, u_xz]    = evaluateSecondOrderDerivativesImpl(this->l(), lfArgs); // (18)!
          const auto &[v_xy, v_xz]    = evaluateSecondOrderDerivativesImpl(this->r(), lfArgs);

          const auto alonguArgs             = replaceAlong(lfArgs, along(u));
          const auto alongvArgs             = replaceAlong(lfArgs, along(v));
          const auto argsForDyzalongv_xArgs = replaceAlong(argsForDyz, along(v_x)); // (19)!
          const auto argsForDyzalongu_xArgs = replaceAlong(argsForDyz, along(u_x));

          const auto u_xyzAlongv = evaluateDerivativeImpl(this->l(), alongvArgs); // (20)!
          const auto v_xyzAlongu = evaluateDerivativeImpl(this->r(), alonguArgs);
          const auto u_yzAlongvx = evaluateDerivativeImpl(this->l(), argsForDyzalongv_xArgs); // (21)!
          const auto v_yzAlongux = evaluateDerivativeImpl(this->r(), argsForDyzalongu_xArgs);

          return Ikarus::eval(u_xyzAlongv + transpose(u_xy) * v_z + transpose(u_xz) * v_y + v_yzAlongux + u_yzAlongvx
                              + transpose(v_xz) * u_y + transpose(v_xy) * u_z + v_xyzAlongu);
        } else if constexpr (LFArgs::hasOneSpatialAll) { // (22)!
          const auto &alongMatrix = std::get<0>(lfArgs.alongArgs.args); // (23)!

          const auto uTimesA = eval(u * alongMatrix);
          const auto vTimesA = eval(v * alongMatrix);

          const auto &[gradu, u_c0, u_c1]  = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs); // (24)!
          const auto &[gradv, v_c0, v_c1]  = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs);
          const auto &[gradu_c0, gradu_c1] = evaluateSecondOrderDerivativesImpl(this->l(), lfArgs); // (25)!
          const auto &[gradv_c0, gradv_c1] = evaluateSecondOrderDerivativesImpl(this->r(), lfArgs);

          const auto graduTimesA = (gradu * alongMatrix.transpose()).eval();
          const auto gradvTimesA = (gradv * alongMatrix.transpose()).eval();

          const auto argsForDyz = lfArgs.extractSecondWrtArgOrFirstNonSpatial();

          const auto alonguAArgs          = replaceAlong(lfArgs, along(uTimesA));
          const auto alongvAArgs          = replaceAlong(lfArgs, along(vTimesA));
          const auto alonggraduTimesAArgs = replaceAlong(argsForDyz, along(graduTimesA));
          const auto alonggradvTimesAArgs = replaceAlong(argsForDyz, along(gradvTimesA));

          const auto u_xyzAlongv            = evaluateDerivativeImpl(this->l(), alongvAArgs);
          const auto v_xyzAlongu            = evaluateDerivativeImpl(this->r(), alonguAArgs);
          const auto v_c0c1AlongGraduTimesA = evaluateDerivativeImpl(this->r(), alonggraduTimesAArgs);
          const auto u_c0c1AlongGradvTimesA = evaluateDerivativeImpl(this->l(), alonggradvTimesAArgs);
          decltype(eval(u_xyzAlongv)) res;

          res = u_xyzAlongv + v_xyzAlongu + v_c0c1AlongGraduTimesA + u_c0c1AlongGradvTimesA;
          for (int i = 0; i < gridDim; ++i)
            res += (transpose(u_c1) * gradv_c0[i] + transpose(v_c1) * gradu_c0[i] + transpose(v_c0) * gradu_c1[i]
                    + transpose(u_c0) * gradv_c1[i])
                   * alongMatrix(0, i);

          return res;

        }
      }
    }
```

1. Compile-time branch for first-order derivatives
2. Evaluates the derivative of `this->l()` in relation to the only derivative contained within the local function arguments `lfArgs`.
3. Evaluates the return value; derivatives and function values are combined as dictated by the product rule.
4. Compile-time branch for second-order derivatives
5. We have four cases for evaluating functions because we are in the second-order derivatives branch: the function value, the function
   derivative w.r.t. the first argument or the second argument, and the function's derivative w.r.t. both arguments.
   Here, the function `evaluateFirstOrderDerivativesImpl` returns the derivatives w.r.t. the first argument and the second argument.
   If we consider the left function as $\boldsymbol{u}(\boldsymbol{\xi},\boldsymbol{u}_I)$, this calls returns
   \begin{flalign*}
   \verb+u_x+ &= \frac{\partial\boldsymbol{u}}{\partial\boldsymbol{\xi}} \quad \text{or} \quad \verb+u_x+ = \frac{\partial\boldsymbol{u}}
   {\partial\xi_0} \quad \text{or} \quad \verb+u_x+ = \frac{\partial\boldsymbol{u}}{\partial\boldsymbol{u}_I}\\\\
   \verb+u_y+ &= \frac{\partial\boldsymbol{u}}{\partial\boldsymbol{u}_J}
   \end{flalign*}
   The first one would be returned if it is called as

   ```cpp
       u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i)));
   ```

   and the second one if it is called as

   ```cpp
       u.evaluateDerivative(gpIndex, wrt(spatial(0),coeff(i)));
   ```

   and the third without any spatial derivative is returned using

   ```cpp
       u.evaluateDerivative(gpIndex, wrt(coeff(i,j)));
   ```

   Therefore, this function separates the two `wrt` arguments and returns the corresponding first order derivatives.
6. Compile-time branch for the case where no spatial derivatives are requested but only the derivatives w.r.t. coefficients are needed.
7. Creates a new argument variable where the `along` argument is replaced by `v`.
8. This function evaluates the derivatives of `l` w.r.t. to both `wrt` arguments. Furthermore, it takes the `along` argument since
   otherwise the returned object would be a 3-dimensional array.
   If we consider the left function as $\boldsymbol{u}(\boldsymbol{\xi},\boldsymbol{u}_I)$  and $\boldsymbol{v}$ of the same size as
   $\boldsymbol{u}$ this calls returns
   \begin{flalign*}
   \verb+u_xyAlongv + &= \frac{\partial^2 u_i }{\partial\boldsymbol{u}_I\partial\boldsymbol{u}_J} v_i
   \end{flalign*}
   This is the same if the following is called:

   ```cpp
       u.evaluateDerivative(gpIndex, wrt(coeff(i,j)),along(v));
   ```

9. Compile-time branch for the case where one spatial derivative and one derivative w.r.t. the coefficients is needed.
10. This function evaluates the derivatives of `l` w.r.t. to both `wrt` arguments.
   If we consider the left function as $\boldsymbol{u}(\boldsymbol{\xi},\boldsymbol{u}_I)$  this calls returns
   \begin{flalign*}
   \verb+u_xy+ &= \frac{\partial^2 \boldsymbol{u} }{\partial\boldsymbol{\xi}\partial\boldsymbol{u}_I}
   \quad \text{or} \quad \verb+u_xy+ = \frac{\partial^2 \boldsymbol{u} }{\partial\xi_0\partial\boldsymbol{u}_I}
   \end{flalign*}
   The first one would be returned if it is called as

    ```cpp
       u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i)));
    ```

    and the second one if it is called as

    ```cpp
    u.evaluateDerivative(gpIndex, wrt(spatial(0),coeff(i)));
    ```

    In the first case the result is stored in an array. Thus, in the first index, the derivative w.r.t. to the first spatial coordinate is stored.
    Therefore, we have

    ```cpp
    spatialAllCoeffDeriv = u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i)));
    spatialAllCoeffDeriv[0] // derivative as in u.evaluateDerivative(gpIndex, wrt(spatial(0),coeff(i)));
    spatialAllCoeffDeriv[1] // derivative as in u.evaluateDerivative(gpIndex, wrt(spatial(1),coeff(i)));
    ```

11. Compile-time branch for the case where one single spatial derivatives and one derivative w.r.t. coefficients is needed.
12. Compile-time branch for the case where all spatial derivatives and one derivative w.r.t. coefficients is needed.
13. The return type here is an array of single spatial derivatives and each derived w.r.t. the coefficient. Thus, the type inside the
    array must be deduced here.
14. Compile-time branch for third order derivatives
15. Compile-time branch for single spatial derivatives
16. To obtain derivatives w.r.t. to the second and third `wrt` argument, we extract the arguments here. E.g., for the following request:

    ```cpp
    u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i,j)),along(matrix));
    ```

    This call would extract the arguments as

    ```cpp
    newArgs =  "wrt(coeff(i,j)),along(matrix))" //THIS IS NO VALID SYNTAX
    ```

    This can then be used as

    ```cpp
    u.evaluateDerivative(gpIndex, newArgs);
    ```

17. As in the second order derivative case, it returns all the three first order derivatives. E.g., for the case,

    ```cpp
    u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i,j)),along(matrix));
    ```

    the returned values would be
    \begin{flalign*}
    \verb+u_x+ &= \frac{\partial \boldsymbol{u} }{\partial\boldsymbol{\xi}} \\\\
    \verb+u_y+ &= \frac{\partial \boldsymbol{u} }{\partial\boldsymbol{u}_I} \\\\
    \verb+u_z+ &= \frac{\partial \boldsymbol{u} }{\partial\boldsymbol{u}_J}
    \end{flalign*}
18. This returns the derivatives w.r.t. the given spatial direction and w.r.t. the first and second coefficient.  E.g., for the case,

    ```cpp
    u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i,j)),along(matrix));
    ```

    the returned values would be
    \begin{flalign*}
    \verb+u_xy+ &= \frac{\partial^2 \boldsymbol{u} }{\partial\boldsymbol{\xi}\partial\boldsymbol{u}_I} \\\\
    \verb+u_xz+ &= \frac{\partial^2 \boldsymbol{u} }{\partial\boldsymbol{\xi}\partial\boldsymbol{u}_J} \\\\
    \end{flalign*}
19. Creates a new argument variable where the `along` argument is replaced by `v_x`.
20. This return call would be

    ```cpp
    u.evaluateDerivative(gpIndex, wrt(spatial(0),coeff(i,j),along(v));
    ```

    In mathematical notation, it returns
    \begin{flalign*}
    \verb+u_xyzAlongv  + &= \frac{\partial^3 u_i }{\partial \xi_0\partial\boldsymbol{u}_I\partial\boldsymbol{u}_J} v_i
    \end{flalign*}
21. This return would be

    ```cpp
    v_x = v.evaluateDerivative(gpIndex, wrt(spatial(0));
    u_yzAlongvx = u.evaluateDerivative(gpIndex, wrt(coeff(i,j),along(v_x));
    ```

    In mathematical notation, it returns
    \begin{flalign*}
    \verb+u_yzAlongvx+ &= \frac{\partial^2 u_i }{\partial\boldsymbol{u}_I\partial\boldsymbol{u}_J} \left[\frac{\partial \boldsymbol{v}}{\xi_0}\right]_i
    \end{flalign*}
22. Compile-time branch for all spatial derivatives
23. Obtain the `along` argument defined, for example, in

    ```cpp
    u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i,j),along(matrix));
    ```

24. Similar to the single spatial case
25. Similar to the single spatial case

If your expression is working, it can be added to `dune/localfefunctions/expressions.hh` by submitting a PR
to [dune-localfefunctions](https://github.com/ikarus-project/dune-localfefunctions).

\bibliography
