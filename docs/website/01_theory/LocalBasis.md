# Local Basis

Each finite element does have some kind of local basis in terms of ansatz functions.
These ansatz function need to be evaluated at the parameter domain of the finite element.

## Interface
Local basis provide the following interface
```cpp
void evaluateFunction(const DomainType& local, Eigen::VectorX<RangeFieldType>& N);
void evaluateDerivative(const DomainType& local,Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& dN);
const Eigen::VectorX<RangeFieldType>& evaluateFunction(const unsigned int& integrationPointIndex);
const Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& evaluateDerivative(const unsigned int& integrationPointIndex);
auto viewOverIntegrationPoints(); // (1)

template <typename IntegrationRule, typename... Ints>
void bind(IntegrationRule&& p_rule, Derivatives<Ints...>&& ints);
```

1. This return a vector of structs of the integration point and its index. Therefore the syntax is usually `#!cpp for (const auto& [gpIndex, gp] : localFunction.viewOverIntegrationPoints()) {...}`

The first two function calls of `evaluateFunction`  and `evaluateDerivative` can be used to calculate the function values 
$ N(\boldsymbol{\xi}) $ and the spatial derivatives $ N_{,\boldsymbol{\xi}}(\boldsymbol{\xi}) $
