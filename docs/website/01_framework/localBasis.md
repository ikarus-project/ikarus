# Local Basis

Each finite element includes a kind of local basis in terms of ansatz functions.
These ansatz functions need to be evaluated in the parameter domain of the finite element.

## Interface

In the exported module ```dune-localfefunctions```, we provide  a thin wrapper with basic caching functionality for dune bases.
For the definitions of bases, refer Chapter 8.2.1 of the Dune book[@sander2020dune].
Thus ```Dune::CachedLocalBasis``` (`#include <dune/cachedlocalBasis/cachedlocalBasis.hh>`) provides the following interface:

   ```cpp
   Dune::LocalBasis(const DuneLocalBasis& p_basis);// (1)!
   void evaluateFunction(const DomainType& local, Eigen::VectorX<RangeFieldType>& N);
   void evaluateJacobian(const DomainType& local,Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& dN);
   void evaluateFunctionAndJacobian(const DomainType& local,Eigen::VectorX<RangeFieldType>& N,
                                    Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& dN);
   
   void evaluateSecondDerivatives(const DomainType& local, Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim*(gridDim + 1) / 2>& dN);
   
   const Eigen::VectorX<RangeFieldType>& evaluateFunction(const unsigned int& integrationPointIndex);
   const Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& evaluateJacobian(const unsigned int& integrationPointIndex);
   const Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim*(gridDim + 1) / 2>& evaluateSecondDerivatives(const unsigned int& integrationPointIndex);
   auto viewOverIntegrationPoints();// (2)!
   auto viewOverFunctionAndJacobian();// (3)!
   
   template <typename IntegrationRule, typename... Ints>
   void bind(IntegrationRule&& p_rule, Derivatives<Ints...>&& ints);
   
   bool isBound(int i) const;// (4)!
   const Dune::QuadraturePoint<DomainFieldType, gridDim>& indexToIntegrationPoint(int i) const;// (5)!
   ```

1. The constructor only accepts a local basis that satisfies the concept `Concepts::DuneLocalBasis`. In keeping with the spirit of
   duck-typing, this also allows for the use of a local basis from the Dune module.
2. This returns a vector of structs containing the integration point and its index. Therefore, the syntax is usually `#!cpp for (const
   auto& [gpIndex, gp] : localFunction.viewOverIntegrationPoints()) {...}`
3. Returns a view over the  ansatz functions and the ansatz function Jacobians at the integration points
4. Checks if the i-th derivatives are bounded.
5. Convert an integration point index to a full integration point.

The first two function calls, `evaluateFunction`  and `evaluateJacobian`, can be used to calculate the function values
\( N(\boldsymbol{\xi}) \) and the spatial derivatives \( N_{,\boldsymbol{\xi}}(\boldsymbol{\xi}) \). One must allocate the objects and
pass them as mutable references.
The same holds for `evaluateSecondDerivatives`.

In contrast to this, there are three other methods that receive an integration point index.
These methods return a `const` reference to the evaluated ansatz function values and their first and second derivatives.

This functionality is dependent on a previous call to `bind(...)`. This binds the local basis to one quadrature rule and caches the passed `bindDerivatives(..)`.
An error is thrown if `#!cpp evaluateFunction(const unsigned int& integrationPointIndex)` is called before binding.
Finally, to bind to an integration rule and cache the value and ansatz function Jacobian, one would use the following syntax:

=== "Usage with integration point index"

    ```cpp linenums="1"
    using LocalFEType = LagrangeSimplexLocalFiniteElement<double,double,dim,2>;
    LocalFEType localFE;

    Dune::CachedLocalBasis(localFE.localBasis());
    const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView.element().type(), order);
    localBasis.bind(rule, bindDerivatives(0, 1));

    for (const auto& [gpIndex, gp] : localBasis.viewOverIntegrationPoints()) {
      const auto& N = localBasis.evaluateFunction(gpIndex);
      const auto& dN = localBasis.evaluateJacobian(gpIndex);
    }
    ```

=== "using integration point coordinates"

    ```cpp linenums="1"
    using LocalFEType = LagrangeSimplexLocalFiniteElement<double,double,dim,2>;
    LocalFEType localFE;

    Dune::CachedLocalBasis(localFE.localBasis());
    const auto& rule = Dune::QuadratureRules<double, Traits::mydim>::rule(localView.element().type(), order);
    Eigen::VectorXd N;
    Eigen::Matrix<double, Eigen::Dynamic, gridDim> dN;

    for(auto& gp : rule){
      localFunction.evaluateFunction(gp.position(), N);
      localFunction.evaluateJacobian(gp.position(), dN);
      localFunction.evaluateFunctionAndJacobian(gp.position(), N, dN);// (1)!
    }
    ```

    1. Alternative to the two lines above (Line 6 and 7)
