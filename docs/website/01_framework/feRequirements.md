# FE requirements

Finite element requirements are a simple way to communicate the needs and expectations of a finite element.

FE requirements are used to pass information from [assemblers](assembler.md) to finite elements.

## Construction and usage

Usually the construction is as follows:

```cpp linenums="1"
FErequirements req = fe.createRequirement()
                           .insertGlobalSolution( d)
                           .insertParameter( lambda);
MatrixType A = sparseFlatAssembler.matrix(req,MatrixAffordance::stiffness,DBCOption::Full);
```

All the methods return a reference to `FErequirements`, so they can be chained together.

As in line 2, the finite element solution can also be inserted. The solution is passed with the `enum` type
`FESolutions::displacement` vector `d`. This stores a reference to the vector.

Additionally, if some parameters are to be passed, use the method `insertParameter` (line 3), where, similar to the
global solutions, an `enum` type of `FEParameter::loadfactor` is passed to indicate the meaning of the parameter,
followed by its value.

Finally, the method `addAffordance` is used to indicate the request required from the finite element.
Thus, there exist scalar, vector, matrix and general affordance collections.

Currently, the following are defined:

{{ inputcpp('ikarus/finiteelements/ferequirements.hh',False,17,52) }}

Inside the finite element, the information can then be conveniently extracted:

```cpp linenums="1"
const auto& d      = req.globalSolution();
const auto& lambda = req.parameter();
```

Thus, the local finite element can be developed.

!!! hint "Affordance"
        It is good style to indicate that you cannot fulfill an affordance by throwing an appropriate exception!

### FE result requirements

!!! bug "Bug"

    This section is outdated!

In addition to the above-mentioned finite element requirements, there are also result requirements.
They accept the same parameter types as the `FErequirements` and add one more, the `ResultType`.
These are used in the `calculateAt` method of [finite elements](finiteElements.md).
They are used to communicate the results required from the finite elements.

Its construction is shown below:

```cpp linenums="1"
ResultRequirements resultRequirements = Ikarus::ResultRequirements()
        .insertGlobalSolution(FESolutions::displacement, d)
        .insertParameter(FEParameter::loadfactor, lambda)
        .addResultRequest(ResultType::cauchyStress,,ResultType::director);
```

The current supported results are:

{{ inputcpp('ikarus/finiteelements/ferequirements.hh',False,53,64) }}

The interface for result requirements is similar to the finite element requirements.
They do, however, support querying specific results to be calculated.

```cpp
if( req.isResultRequested( ResultType::cauchyStress)) {
  ...
}

if( req.isResultRequested( ResultType::BField)) {
  ...
}

if( req.isResultRequested( ResultType::director)) {
  ...
}
```

\bibliography
