<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

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
template<std::size_t ID=0> 
constexpr int order(Dune::index_constant<ID> ); // (2) 
template<std::size_t ID=0> 
auto basis(Dune::index_constant<ID> ); // (3) 
template<std::size_t ID=0> 
auto coefficientsRef(Dune::index_constant<ID>); // (4) 
 
template <typename IntegrationRule, typename... Ints> 
void bind(IntegrationRule&& p_rule, Derivatives<Ints...>&& ints); // (5) 
 
auto clone (); // (6) 
template<typename ScalarType, std::size_t ID=0> 
auto rebindClone (ScalarType, Dune::index_constant<ID>); // (7) 
``` 

1. This returns a vector of structs of the integration point and its index. Therefore the syntax is usually `#!cpp for (const auto& [gpIndex, gp] : localFunction.viewOverIntegrationPoints()) {...}`
2. Return the order of the local function wrt. the coefficients. An id tag can be passed which returns the order wrt a tagged function. For details see [Tagging leaf local functions](#tagging-leaf-local-functions).
3. Return the basis of the local function. An id tag can be passed which returns the basis of a specific tagged function. For details see [Tagging leaf local functions](#tagging-leaf-local-functions).
4. Returns a reference to the coefficient of the underlying leaf local finite elements. An id tag can be passed which returns the basis of a specific tagged function. It can return const and non-const reference. The non-const version is deactivated, if there are more than one leaf node with the passed id tag.  For details see [Tagging leaf local functions](#tagging-leaf-local-functions).
5. This function is passed through to the given `localBasis`. See [Link](localBasis.md)
6. Clones the local function and stores a copy of all leave nodes.
7. Clones the local function and rebinds the scalar type of the coefficients with id tag ID. This becomes hand, if you want to replace doubles with an autodiff type.

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

If we want to mimize this energy w.r.t. the coefficients of the nodes, we need to calculate the energy, gradient and the Hessia w.r.t. the coefficients.  
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
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)), transformWith(Jinv)); 
      Eigen::Vector<double, size> tmp; 
      tmp.setZero(); 
      for (int k = 0; k < gridDimension; ++k) 
        tmp += graduDCoeffs[k] * gradu.col(k);  // (1) 
      g.segment<size>(i * size) += tmp * ("weight from integration point and geo.integrationElement"); 
    } 
  } 
} 
``` 

1. `graduDCoeffs` contains in `graduDCoeffs[0]` the derivatives w.r.t.the coefficient of the first column and at `[1]` w.r.t.the second column of `gradu`

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
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)), transformWith(Jinv)); 
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
In the following table $N^i(\boldsymbol{\xi})$ are the ansatz functions.

| Name                      | Interpolation formula                                         | Note                                                                                                                                                                                                                                                      | Header |
|:--------------------------|:--------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--|
| Standard                    | $$ \boldsymbol{x} = \sum_{i=1}^n N^i(\boldsymbol{\xi}) \boldsymbol{x}_i  $$     |                                                                                                                                                                                                                                                           | `standardLocalFunction.hh`|
| Projection-Based[@grohs_ProjectionBasedFinite2019] | $$ \boldsymbol{x} = P\left(\sum_{i=1}^n N^i(\boldsymbol{\xi}) \boldsymbol{x}_i \right) $$ | This is one version of geometric finite elements. These are finite elements suited for interpolation on manifolds. Here $P: \mathbb{R}^m \rightarrow \mathcal{M}$ is an operator that projects <br /> the usual linear interpolation onto some manifold | `projectionBasedLocalFunction.hh`|

## How to implement your own local functions
If you are interested in implementing your own local function we have prepared the file
[`ikarus/localFunctions/impl/localFunctionTemplate.hh`](https://github.com/IkarusRepo/Ikarus/src/include/ikarus/localFunctions/impl/localFunctionTemplate.hh).

You can copy the file rename the class to your preferred name and then implement the following functions. If you don't need a function you need to delete the corresponding function.
Then if someone calls the corresponding derivative returns a zero matrix.

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
 
 
CoeffDerivMatrix  
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
7. This is called by `localFunction.evaluateDerivative(..., wrt(spatialAll,coeff(j,k)), along(A))`.
8. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i),coeff(j,k)), along(v))`.

## Expressions
We use expression templates[^et] to combine existing local functions to obtain new nested ones.
[^et]:  [Expression templates](https://de.wikipedia.org/wiki/Expression_Templates) are usually used in linear algebra libraries, e.g. [Eigen](https://eigen.tuxfamily.org) or [Blaze](https://bitbucket.org/blaze-lib/blaze/src/master/).
The syntax is similar to the one provided by [UML](https://fenics.readthedocs.io/projects/ufl/en/latest/manual/form_language.html) but only acts on local functions.

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

Currently, we support binary and unary expressions. The following expressions are defined:

| Name        | Mathematical formula                                | Code                        | Note                                                                                  |  
|:------------|:----------------------------------------------------|:----------------------------|:--------------------------------------------------------------------------------------| 
| Sum         | $$ \boldsymbol{f} + \boldsymbol{g}  $$              | `#!cpp f+g`                 | $\boldsymbol{f}$  and  $\boldsymbol{g}$ need to be the same size.                     | 
| DotProduct  | $$ \boldsymbol{f} \cdot \boldsymbol{g} = f_i g_i $$ | `#!cpp dot(f,g)`            | $\boldsymbol{f}$  and  $\boldsymbol{g}$ need to be the same size.                     | 
| normSquared | $$ \boldsymbol{f} \cdot \boldsymbol{f} = f_i f_i $$ | `#!cpp normSquared(f)`      |                                                                                       | 
| Negate      | $$ -\boldsymbol{f}  $$                              | `#!cpp -f`                  |                                                                                       | 
| sqrt        | $$ \sqrt{f}  $$                                     | `#!cpp sqrt(f)`             | The function $f$ needs a scalar return type.                                          | 
| log         | $$ \log{f}  $$                                      | `#!log log(f)`              | The function $f$ needs a scalar return type. Log is the natural logarithm.            | 
| pow         | $$ f^n  $$                                          | `#!cpp pow<n>(f)`           | The function $f$ needs a scalar return type. $n$ is an integer given at compile time. | 
| Scale       | $$  a f , \quad a \in  \mathbf{R}$$                 | `#!cpp a*f` and `#!cpp f/a` | `#!cpp a` has to satisfy `#!cpp std::is_arithmetic<..>`                               | 

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
2. Returns the number of children. 2 for binary expressions and 1 for unary expressions.

!!! note
    To use these expression you can simply include the header by `#!cpp #include <ikarus/localFunctions/expressions.hh>`.

# Tagging leaf local functions
In the context of mixed finite elements. There are usually several local functions that contribute to the energy. These steems from different local basis.
For example consider the Q1P0 element where displacements are interpolated by using the four bilinear ansatz function and the the element-wise constant pressure field.

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


## Writing your own expression
You can also write your own expressions. For this you can look into existing expressions. Especially the sqrt expression and the normSquared expression are the most general unary and binary expression

# Implementing the return value
If you want to implement your own expression you first have to implement the return value.
This is done using the function
```cpp 
template <typename LFArgs> 
auto evaluateValueOfExpression(const LFArgs &lfArgs) const; 
``` 
!!! warning
The interface dictates that the return value needs to be an Eigen type. Thus, even if you want to return a scalar `#!cpp double` you have to wrap it in `#!cpp Eigen::Vector<double, 1>`

Additionally you also have to implement the derivative evaluation. This is done by implementing
```cpp 
template <int DerivativeOrder, typename LFArgs> 
auto evaluateDerivativeOfExpression(const LFArgs &lfArgs) const; 
``` 

# Evaluate underlying functions
Expression always act on already given expression. Therefore, to return the correct quantity for your expression you have to evaluate the underlying quantities.
If you have a unary function you have access to expression using  `#!cpp this->m()` and for binary expressions this is `#!cpp this->l()` and `#!cpp this->r()`.

To evaluate these functions you can use the following syntax.

```cpp 
const auto mEvaluated = evaluateFunctionImpl(this->m(), lfArgs); // (1) 
``` 

1. The syntax is the same for binary expression, e.g.
   ```cpp 
   const auto l_Evaluated = evaluateFunctionImpl(this->l(), lfArgs); 
   const auto r_Evaluated = evaluateFunctionImpl(this->r(), lfArgs); 
   ``` 

The expression fulfill the syntax of a local function thus also derivative can be evaluated.

In the function `evaluateDerivativeOfExpression` the derivative order that the user wants is encoded in the template argument `DerivativeOrder`.
Additionally, the derivative types can also accessed using the booleans

```cpp 
    static constexpr bool hasTwoCoeff; 
    static constexpr bool hasSingleCoeff; 
    static constexpr bool hasNoCoeff; 
    static constexpr bool hasNoSpatial; 
    static constexpr bool hasOneSpatialAll; 
    static constexpr bool hasOneSpatialSingle; 
    static constexpr bool hasOneSpatial; 
``` 

Using the dotproduct as binary example expression we have

```cpp 
template <int DerivativeOrder, typename LFArgs> 
    auto evaluateDerivativeOfExpression(const LFArgs &lfArgs) const { 
      const auto u = evaluateFunctionImpl(this->l(), lfArgs); 
      const auto v = evaluateFunctionImpl(this->r(), lfArgs); 
      if constexpr (DerivativeOrder == 1)  // (1) 
      { 
        const auto u_x = evaluateDerivativeImpl(this->l(), lfArgs); // (2) 
        const auto v_x = evaluateDerivativeImpl(this->r(), lfArgs); 
        return Ikarus::eval(v.transpose() * u_x + u.transpose() * v_x); // (3) 
      } else if constexpr (DerivativeOrder == 2) {   // (4) 
        const auto &[u_x, u_y] = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs); // (5) 
        const auto &[v_x, v_y] = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs); 
        if constexpr (LFArgs::hasNoSpatial and LFArgs::hasTwoCoeff) { // (6) 
          const auto alonguArgs = replaceAlong(lfArgs, along(v)); // (7) 
          const auto alongvArgs = replaceAlong(lfArgs, along(u));  
 
          const auto u_xyAlongv = evaluateDerivativeImpl(this->l(), alongvArgs); // (8) 
          const auto v_xyAlongu = evaluateDerivativeImpl(this->r(), alonguArgs); 
 
          return Ikarus::eval(u_xyAlongv + transpose(u_x) * v_y + transpose(v_x) * u_y + v_xyAlongu); 
        } else if constexpr (LFArgs::hasOneSpatial and LFArgs::hasSingleCoeff) { // (9) 
          const auto u_xy = evaluateDerivativeImpl(this->l(), lfArgs); // (10) 
          const auto v_xy = evaluateDerivativeImpl(this->r(), lfArgs); 
          if constexpr (LFArgs::hasOneSpatialSingle and LFArgs::hasSingleCoeff) { // (11) 
            return Ikarus::eval(transpose(v) * u_xy + transpose(u_x) * v_y + transpose(v_x) * u_y 
                                + transpose(u) * v_xy); 
          } else if constexpr (LFArgs::hasOneSpatialAll and LFArgs::hasSingleCoeff) { // (12) 
            std::array<std::remove_cvref_t<decltype(Ikarus::eval(transpose(v) * u_xy[0]))>, gridDim> res; // (13) 
            for (int i = 0; i < gridDim; ++i) 
              res[i] = Ikarus::eval(transpose(v) * u_xy[i] + transpose(u_x.col(i)) * v_y + transpose(v_x.col(i)) * u_y 
                                    + transpose(u) * v_xy[i]); 
            return res; 
          } 
        } 
      } else if constexpr (DerivativeOrder == 3) { // (14)                                      
        if constexpr (LFArgs::hasOneSpatialSingle) {  // (15) 
          const auto argsForDyz = lfArgs.extractSecondWrtArgOrFirstNonSpatial(); // (16) 
 
          const auto &[u_x, u_y, u_z] = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs); // (17) 
          const auto &[v_x, v_y, v_z] = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs); 
          const auto &[u_xy, u_xz]    = evaluateSecondOrderDerivativesImpl(this->l(), lfArgs); // (18) 
          const auto &[v_xy, v_xz]    = evaluateSecondOrderDerivativesImpl(this->r(), lfArgs); 
 
          const auto alonguArgs             = replaceAlong(lfArgs, along(u)); 
          const auto alongvArgs             = replaceAlong(lfArgs, along(v)); 
          const auto argsForDyzalongv_xArgs = replaceAlong(argsForDyz, along(v_x)); // (19) 
          const auto argsForDyzalongu_xArgs = replaceAlong(argsForDyz, along(u_x)); 
 
          const auto u_xyzAlongv = evaluateDerivativeImpl(this->l(), alongvArgs); // (20) 
          const auto v_xyzAlongu = evaluateDerivativeImpl(this->r(), alonguArgs); 
          const auto u_yzAlongvx = evaluateDerivativeImpl(this->l(), argsForDyzalongv_xArgs); // (21) 
          const auto v_yzAlongux = evaluateDerivativeImpl(this->r(), argsForDyzalongu_xArgs); 
 
          return Ikarus::eval(u_xyzAlongv + transpose(u_xy) * v_z + transpose(u_xz) * v_y + v_yzAlongux + u_yzAlongvx 
                              + transpose(v_xz) * u_y + transpose(v_xy) * u_z + v_xyzAlongu); 
        } else if constexpr (LFArgs::hasOneSpatialAll) { // (22) 
          const auto &alongMatrix = std::get<0>(lfArgs.alongArgs.args); // (23) 
 
          const auto uTimesA = eval(u * alongMatrix); 
          const auto vTimesA = eval(v * alongMatrix); 
 
          const auto &[gradu, u_c0, u_c1]  = evaluateFirstOrderDerivativesImpl(this->l(), lfArgs); // (24) 
          const auto &[gradv, v_c0, v_c1]  = evaluateFirstOrderDerivativesImpl(this->r(), lfArgs); 
          const auto &[gradu_c0, gradu_c1] = evaluateSecondOrderDerivativesImpl(this->l(), lfArgs); // (25) 
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

1. Compile time branch for first order derivatives
2. Evaluates the derivative of the `this->l()` wrt. the only derivative inside the localfunction arguments `lfArgs`.
3. Evaluates the return value and derivatives and function values are combined as dictated by the product rule.
4. Compile time branch for second order derivatives
5. Since we are in the second order derivatives branch, there are 4 case for the evaluation of function.  
   The function value, the function derivative wrt. to the first argument or the second and the function's derivatives wrt. to both arguments.
   Here, the function `evaluateFirstOrderDerivativesImpl` returns the derivatives wrt. to the first argument and the second.
   If we consider the left function as $\boldsymbol{u}(\boldsymbol{\xi},\boldsymbol{u}_I)$ this calls returns
   \begin{flalign*}  
   \verb+u_x+ &= \frac{\partial\boldsymbol{u}}{\partial\boldsymbol{\xi}} \quad \text{or} \quad \verb+u_x+ = \frac{\partial\boldsymbol{u}}{\partial\xi_0} \quad \text{or} \quad \verb+u_x+ = \frac{\partial\boldsymbol{u}}{\partial\boldsymbol{u}_I}\\\\
   \verb+u_y+ &= \frac{\partial\boldsymbol{u}}{\partial\boldsymbol{u}_J}
   \end{flalign*}
   The first one would be returned if the caller uses
   ```cpp 
       u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i))); 
   ``` 
   and the second one
   ```cpp 
       u.evaluateDerivative(gpIndex, wrt(spatial(0),coeff(i))); 
   ``` 
   and the third without any spatial derivative using
   ```cpp 
       u.evaluateDerivative(gpIndex, wrt(coeff(i,j))); 
   ``` 
   Therefore, this function separates the two wrt. arguments and returns the corresponding first order derivatives.
6. Compile time branch for the case where no spatial derivatives are requested bot only wrt. coefficients is needed.
7. Creates a new argument variable where the along argument is replaced by `v`.
8. This function evaluates the derivatives of `l` wrt to both passed wrt arguments. Furthmore, it takes the give along argument since otherwise the returned object would be a 3 dimensional array.
   If we consider the left function as $\boldsymbol{u}(\boldsymbol{\xi},\boldsymbol{u}_I)$  and $\boldsymbol{v}$ of the same size as $\boldsymbol{u}$ this calls returns
   \begin{flalign*}
   \verb+u_xyAlongv + &= \frac{\partial^2 u_i }{\partial\boldsymbol{u}_I\partial\boldsymbol{u}_J} v_i
   \end{flalign*}
   This is the same if the user calls
   ```cpp 
       u.evaluateDerivative(gpIndex, wrt(coeff(i,j)),along(v)); 
   ``` 
9. Compile time branch for the case where one spatial derivatives and one derivative wrt. coefficients is needed.
10. This function evaluates the derivatives of `l` wrt to both passed wrt arguments.
    If we consider the left function as $\boldsymbol{u}(\boldsymbol{\xi},\boldsymbol{u}_I)$  this calls returns
    \begin{flalign*}
    \verb+u_xy+ &= \frac{\partial^2 \boldsymbol{u} }{\partial\boldsymbol{\xi}\partial\boldsymbol{u}_I} \quad \text{or} \quad \verb+u_xy+ = \frac{\partial^2 \boldsymbol{u} }{\partial\xi_0\partial\boldsymbol{u}_I}
    \end{flalign*}
    The first one would be returned if the caller uses
   ```cpp 
       u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i))); 
   ``` 
and the second one
   ```cpp 
   u.evaluateDerivative(gpIndex, wrt(spatial(0),coeff(i))); 
   ``` 
In the first case the result is stored in an array. Thus in the first index the derivative wrt. to the first spatial coordinate is stored.
Therefore we would have in the code
   ```cpp 
   spatialAllCoeffDeriv = u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i))); 
   spatialAllCoeffDeriv[0] // derivative as in u.evaluateDerivative(gpIndex, wrt(spatial(0),coeff(i))); 
   spatialAllCoeffDeriv[1] // derivative as in u.evaluateDerivative(gpIndex, wrt(spatial(1),coeff(i))); 
   ``` 
11. Compile time branch for the case where one single spatial derivatives and one derivative wrt. coefficients is needed.
12. Compile time branch for the case where all spatial derivatives and one derivative wrt. coefficients is needed.
13. The return type here is an array of single spatial derivatives and each derived wrt. the coefficient. Thus the type inside the array must be deduced here.
14. Compile time branch for third order derivatives
15. Compile time branch for single spatial derivatives
16. To obtain derivatives wrt to the second and third wrt argument we extract here the arguments. E.g. if we have the following request
   ```cpp 
    u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i,j)),along(matrix)); 
   ``` 
this call would extract the arguments as
   ```cpp 
    newArgs =  "wrt(coeff(i,j)),along(matrix))" //THIS IS NO VALID SYNTAX 
   ```  
This can be used then as
   ```cpp 
   u.evaluateDerivative(gpIndex, newArgs); 
   ``` 
17. As in the second order derivative case the returns all three first order derivatives. If we would have
    ```cpp 
    u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i,j)),along(matrix)); 
    ``` 
    The returned values would be
    \begin{flalign*}
    \verb+u_x+ &= \frac{\partial \boldsymbol{u} }{\partial\boldsymbol{\xi}} \\\\
    \verb+u_y+ &= \frac{\partial \boldsymbol{u} }{\partial\boldsymbol{u}_I} \\\\
    \verb+u_z+ &= \frac{\partial \boldsymbol{u} }{\partial\boldsymbol{u}_J}  
    \end{flalign*}
18. This returns the derivatives wrt to the given spatial direction and wrt to the first and second coefficient.  If we would have
    ```cpp 
    u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i,j)),along(matrix)); 
    ``` 
    The returned values would be
    \begin{flalign*}
    \verb+u_xy+ &= \frac{\partial^2 \boldsymbol{u} }{\partial\boldsymbol{\xi}\partial\boldsymbol{u}_I} \\\\
    \verb+u_xz+ &= \frac{\partial^2 \boldsymbol{u} }{\partial\boldsymbol{\xi}\partial\boldsymbol{u}_J} \\\\
    \end{flalign*}
19. Creates a new argument variable where the along argument is replaced by `v_x`.
20. This return as the call would be
    ```cpp 
    u.evaluateDerivative(gpIndex, wrt(spatial(0),coeff(i,j),along(v)); 
    ``` 
    In mathematical notation this returns
    \begin{flalign*}
    \verb+u_xyzAlongv  + &= \frac{\partial^3 u_i }{\partial \xi_0\partial\boldsymbol{u}_I\partial\boldsymbol{u}_J} v_i
    \end{flalign*}
22. This return as the call would be
    ```cpp 
    v_x = v.evaluateDerivative(gpIndex, wrt(spatial(0)); 
    u_yzAlongvx = u.evaluateDerivative(gpIndex, wrt(coeff(i,j),along(v_x)); 
    ``` 
    In mathematical notation this returns
    \begin{flalign*}
    \verb+u_yzAlongvx+ &= \frac{\partial^2 u_i }{\partial\boldsymbol{u}_I\partial\boldsymbol{u}_J} \left[\frac{\partial \boldsymbol{v}}{\xi_0}\right]_i
    \end{flalign*}
23. Compile time branch for all spatial derivatives
24. Obtain the along argument give by the caller as in
    ```cpp 
    u.evaluateDerivative(gpIndex, wrt(spatialAll,coeff(i,j),along(matrix)); 
    ``` 
25. As above in the single spatial case
26. As above in the single spatial case

If your expression is working you should add it to `ikarus/localfunctions/expressions.hh`
\bibliography