# Differentiable Function

The class ``DifferentiableFunction`` consists of a collection of a function and its derivatives.

Let us assume a function `f(x)` and its derivative `df(x)`.
Then, a ``DifferentiableFunction`` can be constructed as follows:

```cpp
double x               = 13;
auto fvLambda  = [&](auto&& x) { return f(x); };
auto dfvLambda = [&](auto&& x) { return df(x); };

auto f = Ikarus::makeDifferentiableFunction(functions(fvLambda, dfvLambda), x);
```

!!! note
    It is assumed that the second function is the derivative of the first function, the third function is the derivative of the second
    function (2nd derivative of the first function), and so on.

``functions(...)``  is a helper function. It is necessary to distinguish which argument is a function and which
argument is a parameter.

``f`` provides the following features:

```cpp
void f.operator(x) // (1)!
void derivative(f) // (2)!
```

1. Evaluates the function
2. Returns the function objects corresponding to the first derivative.
