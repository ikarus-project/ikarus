# Non-linear operator

The class ``NonLinearOperator`` provides a collection of a function and its derivatives, including the dependence 
on parameters. 

Let's assume you have a function ``f(x)`` and its derivative ``df(x)``.
Then, a ``NonLinearOperator`` can be constructed as follows:

```cpp
double x               = 13;
auto fvLambda  = [&](auto&& x) { return f(x); };
auto dfvLambda = [&](auto&& x) { return df(x); };

auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(fvLambda, dfvLambda), parameter(x));
```
!!! warning "It is assumed that the second function is the derivative of the first function, the third function is the derivative of the second function (2nd derivative of the first function), etc."

``linearAlgebraFunctions(...)`` and ``parameter(...)`` are helper functions. They are necessary to distinguish which arguments are 
is a function and which argument is parameter. 


``nonLinOp`` now provides the following features:
```cpp
void updateAll() // (1)
void update<n>() // (2)
auto& value() // (3)
auto& derivative() // (4)
auto& secondDerivative() // (5)
auto& nthDerivative<n>() // (6)
auto& lastParameter() // (7)
auto& firstParameter() // (8)
auto& secondParameter() 
auto& nthParameter<n>()
auto subOperator<n,m,...>()

auto nonLinOp = Ikarus::NonLinearOperator(fvLambda, dfvLambda), parameter(x));
```

1. Evaluates all functions
2. Evaluates the n-th function in ``linearAlgebraFunctions(...)`` . Counting starts from 0 as always in C++.
3. Returns the result of the function evaluation.
4. Returns the result of the evaluation of the first derivative (if the function for the first derivative of passed to the nonlinear Operator at construction.)
5. Returns the result of the evaluation of the second derivative (if the function for the second derivative of passed to the nonlinear Operator at construction.).
6. Returns the result of the evaluation of the n-th derivative (if the function for the third derivative of passed to the nonlinear Operator at construction.).
7. 