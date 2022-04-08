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
      localFunction.evaluateDerivative(gpIndex, wrt(spatialall)); // (1)
      localFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv)); // (2)
    }
    ```

    1. Compute the spatial Jacobian of localFunction
    2. Compute the spatial Jacobian of localFunction and transform it to physical coordinates

=== "using integration point coordinates"

    ``` c++
    using namespace Ikarus::DerivativeDirections;
    for(auto& gp : rule){
      localFunction.evaluateDerivative(gp.position(), wrt(spatialall)); // (1)
      localFunction.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv)); // (2)
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

Instead of passing `spatialall` to `wrt(..)`, there are other helper such as

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(spatial(0))); // (1) 
localFunction.evaluateDerivative(gpIndex, wrt(spatial(1))); // (2)
```

1. Compute the first column of the spatial Jacobian of localFunction
2. Compute the second column of the spatial Jacobian of localFunction

which can also be combined with `transformWith(Jinv)`.

## Derivatives w.r.t. coefficients
```cpp
localFunction.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(j));
```
which implements for simple interpolation in vector space valued functions,e.g. $f(\boldsymbol{\xi}) = \sum_{I=1}^n N^I(\boldsymbol{\xi}) \boldsymbol{x}_I$ the following

$$
[\boldsymbol{A}]_{ij}  = A_{ij} =  \frac{\partial f_i(\boldsymbol{\xi})}{\partial \boldsymbol{x}_j}
$$

and the second derivative

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(coeffs,coeffs), along(q), coeffIndices(j,k));
```

$$
[\boldsymbol{B}]_{jk} =  B_{jk} = q_i A_{ijk} =  \frac{\partial^2 (q_i  f_i(\boldsymbol{\xi}))}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}
$$

where $\boldsymbol{q}$ is an arbitrary vector of the same size as $f$, i.e. it is the direction of the derivative in this case. $ \boldsymbol{A} $ and $ \boldsymbol{B} $ is simply the returned matrix and they do not have a special meaning. If we would not pass the vector the result would be a third order tensor for a vector valued function $f$.
Therefore the simply return a matrix. This helps for readablilty and for speed. See the [example](#example-dirichlet-energy) for details.
## Derivatives w.r.t. coefficients and spatial derivatives
Spatial derivatives and derivatives w.r.t. the coefficients can be combined. Therefore, it is legal to call

```cpp
auto B = localFunction.evaluateDerivative(gpIndex, wrt(coeffs,coeffs,spatialall), along(q), coeffIndices(j,k));
auto B0 = localFunction.evaluateDerivative(gpIndex, wrt(coeffs,coeffs,spatial(0)), along(q), coeffIndices(j,k));
auto B1 = localFunction.evaluateDerivative(gpIndex, wrt(coeffs,coeffs,spatial(1)), along(q), coeffIndices(j,k));
```

The first line is then equivalent to

$$
[\boldsymbol{B}]_{ljk} =  B_{ljk} = q_i A_{iljk} =  \frac{\partial^2 ([\operatorname{grad}_\boldsymbol{\xi} f(\boldsymbol{\xi})]_{il} q_i )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
$$

this returns an object where the first index contains the spatial derivative w.r.t. $\xi_0$.
Thus we have 

\begin{align}
\boldsymbol{B}[0]_{jk} = \frac{\partial^2 ([\operatorname{grad}_{\xi^0} f(\xi)]_{i} q_i )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}, \\
\boldsymbol{B}[1]_{jk} = \frac{\partial^2 ([\operatorname{grad}_{\xi^1} f(\xi)]_{i} q_i )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
\end{align}

These objects are also returned when the second and third line above are used.

Again all of these function calls can be combined with `transformWith()` as

```cpp
localFunction.evaluateDerivative(gpIndex, wrt(coeffs,coeffs,spatialall), along(q),transformWith(Jinv), coeffIndices(j,k));
```

which computes

$$
\frac{\partial^2 ([\operatorname{grad}_\boldsymbol{x} f(\boldsymbol{\xi})]_{il} q_i )}{\partial \boldsymbol{x}_j\partial \boldsymbol{x}_k}.
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
    const auto gradu = uFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));
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
    const auto gradu = uFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));
    for (auto i : fe.size()) { //loop over coeffs, i.e.nodes of the finite element
      const auto graduDCoeffs
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(i));
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
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(i));
      for (auto j : fe.size()) { //loop over coeffs, i.e.nodes of the finite element
        const auto graduDCoeffsJ
          = uFunction.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(j));
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
1. This is called by `localFunction.evaluateDerivative(..., wrt(spatialall))`.
2. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i)))`.
3. This is called by `localFunction.evaluateDerivative(..., wrt(coeff))`.
4. This is called by `localFunction.evaluateDerivative(..., wrt(coeff,coeff))`.
5. This is called by `localFunction.evaluateDerivative(..., wrt(spatialall,coeff))`.
6. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i),coeff))`.
7. This is called by `localFunction.evaluateDerivative(..., wrt(spatialall,coeff,coeff))`.
8. This is called by `localFunction.evaluateDerivative(..., wrt(spatial(i),coeff,coeff))`.


\bibliography