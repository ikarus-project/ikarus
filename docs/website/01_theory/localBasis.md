# Local Basis

Each finite element does have some kind of local basis in terms of ansatz functions.
These ansatz function need to be evaluated at the parameter domain of the finite element.

## Interface
Local basis provide the following interface
```cpp
LocalBasis(const DuneLocalBasis& p_basis)  // Constructor (1)
void evaluateFunction(const DomainType& local, Eigen::VectorX<RangeFieldType>& N);
void evaluateJacobian(const DomainType& local,Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& dN);
void evaluateFunctionAndJacobian(const DomainType& local,Eigen::VectorX<RangeFieldType>& N,
                                 Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& dN);

const Eigen::VectorX<RangeFieldType>& evaluateFunction(const unsigned int& integrationPointIndex);
const Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& evaluateJacobian(const unsigned int& integrationPointIndex);
auto viewOverIntegrationPoints(); // (2)

template <typename IntegrationRule, typename... Ints>
void bind(IntegrationRule&& p_rule, Derivatives<Ints...>&& ints);
```

1. Using the concept `Concepts::DuneLocalBasis`  the constructor only accepts local basis that satisfies this concept. This also allows a local basis which behave lik a local basis of dune in the spirit of duck-typing.
2. This return a vector of structs of the integration point and its index. Therefore the syntax is usually `#!cpp for (const auto& [gpIndex, gp] : localFunction.viewOverIntegrationPoints()) {...}`

The first two function calls of `evaluateFunction`  and `evaluateJacobian` can be used to calculate the function values 
\( N(\boldsymbol{\xi}) \) and the spatial derivatives \( N_{,\boldsymbol{\xi}}(\boldsymbol{\xi}) \). The objects where this is stored you have to allocate yourself and have to pass as mutable reference.

In contrast to this there exists two other methods that receive an integration point index. 
These methods return a const reference to the evaluated ansatz function values and derivatives.

This functionality depends on an earlier call to `bind(...)`. This binds the local basis to one quadrature rule and caches the passed `bindDerivatives(..)`.
If one calls `#!cpp evaluateFunction(const unsigned int& integrationPointIndex)` before bind an error is thrown.
Finally, to bind to an integration rule and cache the value and the ansatz function jacobian one would call:

=== "Usage with integration point index"

    ```cpp linenums="1"
    const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
    localBasis.bind(rule, bindDerivatives(0, 1));

    for (const auto& [gpIndex, gp] : localBasis.viewOverIntegrationPoints()) {
      const auto& N = localBasis.evaluateFunction(gpIndex);
      const auto& dN = localBasis.evaluateJacobian(gpIndex);
    }
    ```

=== "using integration point coordinates"

    ```cpp linenums="1"
    const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView_.element().type(), order);
    Eigen::VectorXd N;
    Eigen::Matrix<double, Eigen::Dynamic, gridDim> dN;

    for(auto& gp : rule){
      localFunction.evaluateFunction(gp.position(), N); 
      localFunction.evaluateJacobian(gp.position(), dN);
      localFunction.evaluateFunctionAndJacobian(gp.position(), N, dN); // (1) 
    }
    ```

    1. Alternative to the two lines above (Line 6 and 7)