# Local functions

This section explains the concept of local functions.

Local functions are functions which are bound to single grid elements.
Therefore they are constructed from some local basis and a coefficient vector.

# Theory
Local functions need to be evaluated in the local coordinate system $\mathbb{\xi} \in \mathbb{R}^n$:
$$
f: \boldsymbol{\xi}^n \rightarrow \mathbb{R}^m
$$
## Interface
Local functions provide the following interface
```cpp
  FunctionReturnType evaluateFunction(const DomainType& local);
  FunctionReturnType evaluateFunction(const unsigned int& integrationPointIndex);
  auto evaluateDerivative(const DomainType& local,...);
  auto evaluateDerivative(const unsigned int& integrationPointIndex,...);
  auto viewOverIntegrationPoints(); // (1)
```
1. This return a vector of structs of the integration point and its index. Therefore the syntax is usually `#!cpp for (const auto& [gpIndex, gp] : localFunction.viewOverIntegrationPoints()) {...}`

The $...$ in the `evaluateDerivative` function call are several variadic templates.
In action this looks like
```cpp
  localFunction.evaluateDerivative(gpIndex, wrt(spatialall)); 
  localFunction.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv)); 
```

## Implementations
The following local functions are currently available:

## Implementations
In the following we sumerize the local functions that are currently available.
In the follwing table $N^i(\boldsymbol{\xi})$ are the ansatz functions.

| Name                      | Interpolation formula                                         | Note                                                                                                                                                                                                                                                      | Header |
|:--------------------------|:--------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--|
| Standard                    | $$ \boldsymbol{x} = \sum_{i=1}^n N^i(\boldsymbol{\xi}) \boldsymbol{x}_i $$     |                                                                                                                                                                                                                                                           | `StandardLocalFunction.h`|
| Projection-Based[@grohs_ProjectionBasedFinite2019] | $$ \boldsymbol{x} = P\left(\sum_{i=1}^n N^i(\boldsymbol{\xi}) \boldsymbol{x}_i \right) $$ | This is one version of geometric finite elements. These are finite elements suited for interpolation on manifolds. \ Here $P: \mathbb{R}^m \rightarrow \mathcal{M}$ is an operator that projects <br /> the usual linear interpolation onto some manifold | `ProjectionBasedLocalFunction.h`|


## How to implement your own local functions
If you are interested in implementing your own local function we have prepared the file
`ikarus/LocalFunctions/LocalFunctionTemplate.h`[Link]()

{{ inputcpp('src/include/ikarus/LocalFunctions/LocalFunctionTemplate.h',true) }}

\bibliography