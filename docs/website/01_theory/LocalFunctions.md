# Local functions

This section explains the concept of local functions.

Local functions are functions which are bound to single grid elements.
Therefore they are constructed from some local basis and a coefficient vector.

Usually local functions need to be evaluated in the local coordinate system \( \mathbb{\xi} \in T_{\text{ref}} \subset\mathbb{R}^n \) :

$$
f: \boldsymbol{\xi}^n \rightarrow \mathbb{R}^m
$$

where $T_{\text{ref}}$ is the reference element, e.g. for a cube $T_{\text{ref}}= [0,1]^d$.
## Interface
Local functions provide the following interface
```cpp
FunctionReturnType evaluateFunction(const DomainType& local);
FunctionReturnType evaluateFunction(const unsigned int& integrationPointIndex);
auto evaluateDerivative(const DomainType& local,...);
auto evaluateDerivative(const unsigned int& integrationPointIndex,...);
auto viewOverIntegrationPoints(); // (1)
template<std::size_t ID>
constexpr int order(Dune::index_constant<ID> );
template<std::size_t ID>
auto basis(Dune::index_constant<ID> );



template <typename IntegrationRule, typename... Ints>
void bind(IntegrationRule&& p_rule, Derivatives<Ints...>&& ints); // (2)
```

1. This returns a vector of structs of the integration point and its index. Therefore the syntax is usually `#!cpp for (const auto& [gpIndex, gp] : localFunction.viewOverIntegrationPoints()) {...}`
2. This function is passed through to the given `localBasis`. See [Link](LocalBasis.md)

The "..." in the `evaluateDerivative` function call are several variadic templates.
In action this looks like


=== "Usage with integration point index"

    ``` c++
    using namespace Ikarus::DerivativeDirections;
    localFunction.bind(rule, bindDerivatives(0,1));    
    for(auto& [gpIndex, gp] : localFunction.viewOverIntegrationPoints()){
      localFunction.evaluateDerivative(gpIndex, wrt(spatialAll)); // (1)
      localFunction.evaluateDerivative(gpIndex, wrt(spatialAll), transformWith(Jinv)); // (2)
    }
    ```

    1. Compute the spatial Jacobian of localFunction
    2. Compute the spatial Jacobian of localFunction and transform it to physical coordinates

=== "using integration point coordinates"

    ``` c++
    using namespace Ikarus::DerivativeDirections;
    for(auto& gp : rule){
      localFunction.evaluateDerivative(gp.position(), wrt(spatialAll)); // (1)
      localFunction.evaluateDerivative(gp.position(), wrt(spatialAll), transformWith(Jinv)); // (2)
    }
    ```

    1. Compute the spatial Jacobian of localFunction
    2. Compute the spatial Jacobian of localFunction and transform it to physical coordinates

where the first call implements

$$
\operatorname{grad}_\boldsymbol{\xi} f : \boldsymbol{\xi} \rightarrow \mathbb{R}^{m \times d}.
$$

The second one respect the fact that the local function in reality is defined in some physical space $X$ with the coordinate $\boldsymbol{x}$.
Therefore, it transforms the Jacobian from the reference element $\operatorname{grad}_{\boldsymbol{\xi}}$ to the Jacobian in physical space $\operatorname{grad}_\boldsymbol{x}$ . E.g. it usually implements
 
$$
\operatorname{grad}_\boldsymbol{x} = \operatorname{grad}_{\boldsymbol{\xi}} \boldsymbol{J}^{-1} 
$$

where $J$ is the Jacobian of the mapping from the reference element $T_{\text{ref}}$ to the element living in physical space $T$.
For details see [@sander2020dune] page 22.

Instead of passing `spatialAll` to `wrt(..)`, there are other helper such as

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(spatial(0))); // (1) 
localFunction.evaluateDerivative(gpIndex, wrt(spatial(1))); // (2)
```

1. Compute the first column of the spatial Jacobian of localFunction
2. Compute the second column of the spatial Jacobian of localFunction

which can also be combined with `transformWith(Jinv)`.

## Derivatives w.r.t. coefficients
```cpp
localFunction.evaluateDerivative(gpIndex, wrt(coeff(j)));
```
which implements for a in vector space valued function (steeming from interpolation),e.g. $f(\boldsymbol{\xi}) = \sum_{I=1}^n N^I(\boldsymbol{\xi}) \boldsymbol{x}_I$ the following

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

where $\boldsymbol{q}$ is an arbitrary vector of the same size as $f$, i.e. it is the direction of the derivative in this case. $ \boldsymbol{A} $ and $ \boldsymbol{B} $ is simply the returned matrix and they do not have a special meaning. If we would not pass the vector the result would be a third order tensor for a vector valued function $f$.
Therefore the simply return a matrix. This helps for readablilty and for speed. See the [example](#example-dirichlet-energy) for details.
## Derivatives w.r.t. coefficients and spatial derivatives
Spatial derivatives and derivatives w.r.t. the coefficients can be combined. Therefore, it is legal to call

```cpp
auto B = localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k),spatialAll), along(Q));
auto b1 = localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k),spatial(0)), along(q));
auto b2 = localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k),spatial(1)), along(q));
```

!!! warning 
    The order of spatial and coeff derivatives does not matter. The returned value is always re-arranged that the first derivative is the spatial one.

The first line is then equivalent to

$$
[\boldsymbol{B}]_{jk} =  B_{jk} = Q_{il} A_{iljk} =  \frac{\partial^2 ([\operatorname{grad}_\boldsymbol{\xi} f(\boldsymbol{\xi})]_{il} Q_{il} )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
$$

For the second and third line we have

\begin{align}
\boldsymbol{b}_{0,jk} = \frac{\partial^2 ([\operatorname{grad}_{\xi^0} f(\xi)]_{i} q_i )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}, \\
\boldsymbol{b}_{1,jk} = \frac{\partial^2 ([\operatorname{grad}_{\xi^1} f(\xi)]_{i} q_i )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
\end{align}

These objects are also returned when the second and third line above are used.

Again all of these function calls can be combined with `transformWith()` as

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(coeff(j,k),spatialAll), along(Q),transformWith(Jinv));
```

which computes

$$
\frac{\partial^2 ([\operatorname{grad}_\boldsymbol{x} f(\boldsymbol{\xi})]_{il} Q_{il} )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
$$


!!! warning "Warning"
    Currently only first order spatial derivatives and second order derivatives w.r.t. the coefficients are supported.

## Example Dirichlet energy
This examples shows how the energy, gradient and Hessian of a [dirichlet energy](https://en.wikipedia.org/wiki/Dirichlet_energy) can be calculated.
$$
E(\boldsymbol{u}) = \frac{1}{2} \int_\Omega ||\operatorname{grad}_\boldsymbol{x} \boldsymbol{u}(\boldsymbol{x})|| ^2 \textrm{d} \boldsymbol{x}
$$

If we want to mimize this energy w.r.t. the coefficients of the nodes, we need to calculate the energy, gradient and the Hessia w.r.t. the coefficents. 
Of course this depends on the optimization algorithms, but for now lets keep it simple.

```cpp
auto dirichletEnergy() {
  double energy = 0;
  //... bind localBasis to some integration rule
  // and create uNodalCoeffs
  Ikarus::StandardLocalFunction uFunction(localBasis, uNodalCoeffs);
  for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
    //.. calculate the inverse Jacobian of the geometry
    const auto gradu = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll), transformWith(Jinv));
    energy+= 0.5*(gradu.transpose()*gradu).trace()* ("weight from integration point and geo.integrationElement");
  }
}
```

```cpp
auto gradientDirichletEnergy(Eigen::VectorXd& g) {
  //... bind localBasis to some integration rule
  // and create uNodalCoeffs
  constexpr int size =  // spatial size of u
      Ikarus::StandardLocalFunction uFunction(localBasis, uNodalCoeffs);
  for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
    //.. calculate the inverse Jacobian of the geometry
    const auto gradu = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll), transformWith(Jinv));
    for (auto i : fe.size()) { //loop over coeffs, i.e.nodes of the finite element
      const auto graduDCoeffs
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeffs), transformWith(Jinv), coeffIndices(i));
      Eigen::Vector<double, size> tmp;
      tmp.setZero();
      for (int k = 0; k < gridDimension; ++k)
        tmp += graduDCoeffs[k] * gradu.col(k);  // (1)
      g.segment<size>(i * size) += tmp * ("weight from integration point and geo.integrationElement");
    }
  }
}
```

1. `graduDCoeffs` contains in `graduDCoeffs[0]` the derivatives w.r.t.the coefficient of the first column and at `[1]` w.r.t.the second colum of `gradu`

```cpp
auto hessianDirichletEnergy(Matrix& h) {
  //... bind localBasis to some integration rule
  // and create uNodalCoeffs
  constexpr int size =  // spatial size of u
      Ikarus::StandardLocalFunction uFunction(localBasis, uNodalCoeffs);
  for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
    //.. calculate the inverse Jacobian of the geometry
    for (auto i : loop over coeffs, i.e.nodes of the finite element) {
      const auto graduDCoeffsI
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeffs), transformWith(Jinv), coeffIndices(i));
      for (auto j : fe.size()) { //loop over coeffs, i.e.nodes of the finite element
        const auto graduDCoeffsJ
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeffs), transformWith(Jinv), coeffIndices(j));
        Eigen::Matrix<double, size, size> tmp;
        tmp.setZero();
        for (int k = 0; k < gridDimension; ++k)
          tmp += graduDCoeffsI[k] * graduDCoeffsJ[k];
        h.block<size, size>(i * size, j * size) += tmp * ("weight from integration point and geo.integrationElement");
      }
    }
  }
}
```

## Implementations
In the following we summarize the local functions that are currently available.
In the follwing table $N^i(\boldsymbol{\xi})$ are the ansatz functions.

| Name                      | Interpolation formula                                         | Note                                                                                                                                                                                                                                                      | Header |
|:--------------------------|:--------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--|
| Standard                    | $$ \boldsymbol{x} = \sum_{i=1}^n N^i(\boldsymbol{\xi}) \boldsymbol{x}_i  $$     |                                                                                                                                                                                                                                                           | `standardLocalFunction.hh`|
| Projection-Based[@grohs_ProjectionBasedFinite2019] | $$ \boldsymbol{x} = P\left(\sum_{i=1}^n N^i(\boldsymbol{\xi}) \boldsymbol{x}_i \right) $$ | This is one version of geometric finite elements. These are finite elements suited for interpolation on manifolds. Here $P: \mathbb{R}^m \rightarrow \mathcal{M}$ is an operator that projects <br /> the usual linear interpolation onto some manifold | `projectionBasedLocalFunction.hh`|

## How to implement your own local functions
If you are interested in implementing your own local function we have prepared the file
[`ikarus/LocalFunctions/LocalFunctionTemplate.h`](https://github.com/IkarusRepo/Ikarus/src/include/ikarus/LocalFunctions/LocalFunctionTemplate.h).

You can copy the file rename the class to your preferred name and then implement the following functions. If you don't need a function you need to delete the corresponding function.
Then if someone calls the corresponding derivative the call fails at compile time.

```cpp
FunctionReturnType evaluateEmbeddingFunctionImpl(const Eigen::VectorXd& N) const { return FunctionReturnType{}; } // (0)

Jacobian evaluateDerivativeWRTSpaceAllImpl(const AnsatzFunctionType& N, 
                                           const AnsatzFunctionJacobian& dN) const {...} // (1)

JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const AnsatzFunctionType& N, 
                                                     const AnsatzFunctionJacobian& dN,
                                                     int spaceIndex) const {...} // (2)


CoeffDerivMatrix evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N,
                                                 const AnsatzFunctionJacobian& dN,
                                                 int coeffsIndex) const {...} // (3)

CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffs(const AnsatzFunctionType& N,
                                                   const AnsatzFunctionJacobian&,
                                                   const AlongType& along,
                                                   const std::array<size_t, gridDim>& coeffsIndex) const {...} // (4)

std::array<CoeffDerivMatrix, gridDim> 
        evaluateDerivativeWRTCoeffsANDSpatialImpl(const AnsatzFunctionType& N, 
                                                  const AnsatzFunctionJacobian& dN, 
                                                  int coeffsIndex) const {...} // (5)


CoeffDerivMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                                 const AnsatzFunctionJacobian& dN,
                                                                 const int coeffsIndex, 
                                                                 const int spatialIndex) const {...} // (6)


std::array<CoeffDerivMatrix, gridDim> 
        evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(const AnsatzFunctionType& N, 
                                                               const AnsatzFunctionJacobian& dN, 
                                                               const AlongType& along,
                                                               const std::array<size_t, gridDim>& coeffsIndex
                                                               ) const {...} // (7)
                                                                                             
CoeffDerivMatrix 
        evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(const AnsatzFunctionType& N, 
                                                                     const AnsatzFunctionJacobian& dN, 
                                                                     const AlongType& along,
                                                                     const std::array<size_t, gridDim>& coeffsIndex, c
                                                                     const int spatialIndex) const {...} // (8)
```

0. This is called by `localFunction.evaluateFunction(...)`.
1. This is called by `localFunction.evaluateDerivative(..., wrt(spatialAll))`.
2. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i)))`.
3. This is called by `localFunction.evaluateDerivative(..., wrt(coeff(j)))`.
4. This is called by `localFunction.evaluateDerivative(..., wrt(coeff(j,k)))`.
5. This is called by `localFunction.evaluateDerivative(..., wrt(spatialAll,coeff(j)))`.
6. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i),coeff(j)))`.
7. This is called by `localFunction.evaluateDerivative(..., wrt(spatialAll,coeff(j,k)))`.
8. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i),coeff(j,k)))`.

## Expressions
We use expression templates[^et] to combine existing local functions to obtain new nested ones.
[^et]:  [Expression templates](https://de.wikipedia.org/wiki/Expression_Templates) are usually used in linear algebra libraries, e.g. [Eigen](https://eigen.tuxfamily.org) or [Blaze](https://bitbucket.org/blaze-lib/blaze/src/master/).

For example consider the following code 
```cpp
...
auto f = Ikarus::StandardLocalFunction(localBasis, coeffVectors0);
auto g = Ikarus::StandardLocalFunction(localBasis, coeffVectors1);
```
we create here two local functions that satisfy the interface described above. 
Now it is possible to combine these functions and get an object that also satisfies the concept above.
Thus the following is possible:
```cpp
...
auto f = Ikarus::StandardLocalFunction(localBasis, coeffVectors0);
auto g = Ikarus::StandardLocalFunction(localBasis, coeffVectors1);
auto k = f+g;
k.evaluateDerivative(ipIndex, wrt(coeff(i), spatial(d)));
```

Currently, we support binary and unary expressions. The following expression are defined:

| Name       | Mathematical formula                                | Code                        | Note                                                              | 
|:-----------|:----------------------------------------------------|:----------------------------|:------------------------------------------------------------------|
| Sum        | $$ \boldsymbol{f} + \boldsymbol{g}  $$              | `#!cpp f+g`                 | $\boldsymbol{f}$  and  $\boldsymbol{g}$ need to be the same size. |
| DotProduct | $$ \boldsymbol{f} \cdot \boldsymbol{g} = f_i g_i $$ | `#!cpp dot(f,g)`            |  $\boldsymbol{f}$  and  $\boldsymbol{g}$ need to be the same size.                                                                 |
| Negate     | $$ -\boldsymbol{f}  $$                              | `#!cpp -f`                  |                                                                   |
| sqrt       | $$ \sqrt{f}  $$                                     | `#!cpp sqrt(f)`             | The function $f$ needs a scalar return type.                      |
| Scale      | $$  a f , \quad a \in  \mathbf{R}$$                 | `#!cpp a*f` and `#!cpp f/a` | `#!cpp a` has to satisfy `#!cpp std::is_arithmetic<..>`           |

These expressions can be nested. Thus, it is valid to write something like
```cpp
auto f = Ikarus::StandardLocalFunction(localBasis, coeffVectors0);
auto g = Ikarus::StandardLocalFunction(localBasis, coeffVectors1);
auto k = -sqrt(dot(2*f+f,5*g));
k.evaluateDerivative(ipIndex, wrt(coeff(i), spatial(d)));
```

To use these expression there are addition exported static types for all expressions
```cpp
constexpr bool isLeaf; // (1)
constexpr bool children; // (2)
```

1. This is true if the underlying expression is one of the above Local functions that really contain the coefficients, see [Implementations](#implementations).
2. Returns the number of childs 2 for binary expressions and 1 for unary expressions.

!!! note
    To use these expression you can simply include the header by `#!cpp #include <ikarus/localFunctions/expressions.hh>`.

# Tagging leaf local functions
In the context of mixed finite elements. There are usually several local functions that contribute to the energy. These steems from different local basis.
For example consider the Q1P0 element where displacments are interpolated by using the four bilinear ansatz function and the the element-wise constant pressure field.

Thus we need to differentiate wrt. different coefficients. This can be done by tagging the local function by construction.
```cpp
using namespace Dune::Indices;
auto f = Ikarus::StandardLocalFunction(localBasis0, coeffVectors0,0_);
auto g = Ikarus::StandardLocalFunction(localBasis1, coeffVectors1,1_);
auto k = dot(f,g);
k.evaluateDerivative(ipIndex, wrt(coeff(0_,i,1_,j)));
```
To explain the last line above lets consider that the function f is constructed as $f= \sum_{I=0}^n N^I f_i$ and similar 
$g= \sum_{I=0}^m M^I g_i$, where $N$ and $M$ are some ansatz functions and $f_I$ and $g_I$ are nodal coefficients.

Thus the above call translates to

\begin{align}
\boldsymbol{M}_{0,1}[J,K] = \frac{\partial^2 (f_{i} g_i )}{\partial \boldsymbol{f}_J\partial \boldsymbol{g}_K}.
\end{align}

If we would calculate the complete hessian of $dot(f,g)$ we can do this by

```cpp
using namespace Dune::Indices;
auto hessianDirichletEnergy(Matrix& h) {
  //... bind localBasis to some integration rule
  using namespace Dune::Indices;
  auto f = Ikarus::StandardLocalFunction(localBasis0, coeffVectors0,0_);
  auto g = Ikarus::StandardLocalFunction(localBasis1, coeffVectors1,1_);
  auto k = dot(f,g);
  constexpr int sizef = f.correctionSize; // spatial size of the correction of the coefficients of f
  constexpr int sizeg = g.correctionSize; // spatial size of the correction of the coefficients of g
  constexpr int coeffSizef = coeffVectors0.size();
  constexpr int coeffSizeg = coeffVectors1.size();
  
  Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<MatrixBlock00,MatrixBlock01>,
                                       Dune::MultiTypeBlockVector<MatrixBlock10,MatrixBlock11> > KBlocked; // (1)

  
  for (const auto& [ipIndex, gp] : k.viewOverIntegrationPoints()) {
    for (size_t I = 0; I < coeffSizef; ++I)
      for (size_t J = 0; J < coeffSizef; ++J) 
        KBlocked[0_,0_].block<sizef, sizef>(I * sizef, J * sizef) 
          += k.evaluateDerivative(ipIndex, wrt(coeff(0_,I,0_,J)))* ("weight from integration point and geo.integrationElement");
    
    for (size_t I = 0; I < coeffSizef; ++I)
      for (size_t J = 0; J < coeffSizeg; ++J)
        KBlocked[0_,1_].block<sizef, sizeg>(I * sizef, J * sizeg) 
          += k.evaluateDerivative(ipIndex, wrt(coeff(0_,I,1_,J)))* ("weight from integration point and geo.integrationElement");
    
    for (size_t I = 0; I < coeffSizeg; ++I)
      for (size_t J = 0; J < coeffSizeg; ++J)
        KBlocked[1_,1_].block<sizeg, sizeg>(I * sizeg, J * sizeg) 
          += k.evaluateDerivative(ipIndex, wrt(coeff(1_,I,1_,J)))* ("weight from integration point and geo.integrationElement");
      
    for (size_t I = 0; I < coeffSizeg; ++I)
      for (size_t J = 0; J < coeffSizef; ++J)
        KBlocked[1_,0_].block<sizef, sizeg>(I * sizeg, J * sizef) 
          += k.evaluateDerivative(ipIndex, wrt(coeff(1_,I,0_,J)))* ("weight from integration point and geo.integrationElement");
    }
}
```

1. This Block structure is not necessary. In this example all types (MatrixBlock00,MatrixBlock01,MatrixBlock10,MatrixBlock11) are considered as `#!cpp Eigen::MatrixXd`.

\bibliography