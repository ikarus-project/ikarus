<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de

SPDX-License-Identifier: CC-BY-SA-4.0
-->

# FE requirements

Finite element requirements are a simply way to communicate your needs and expectations from a finite element.

FE requirements are used to pass information from assemblers to finite elements. 
## Construction
Usually the construction is as follows.
```cpp linenums="1"
FErequirements req = FErequirementsBuilder()
                           .insertGlobalSolution(FESolutions::displacement, d)
                           .insertParameter(FEParameter::loadfactor, lambda)
                           .addAffordance(MatrixAffordances::stiffness)
                           .build();
MatrixType A = sparseFlatAssembler.getReducedMatrix(req);
```

As you can see to construct requirements we used the builder pattern[@gamma1995design].

Thus to construct `FErequirements` you create a `FErequirementsBuilder`.
You can then chain your requirements together.

You can insert solution from your finite element solution algorithm as in line 2. There, the type of the soultion is passed with the enum type 
`FESolutions::displacement` with the vector `d`. This stores a reference to the vector.

Additionally, if you have some parameters you want to pass you can call the method `insertParameter` as in line 3, where similar as for the 
global solutions a enum `FEParameter::loadfactor` is passed to indicate the meaning of the parameter and after this the value is passed.

Finally there is the method `addAffordance` which is used to indicate your request what you want from the finite element.
Thus, there are scalar, vector and matrix affordances.

The method `build()` constructs at the end the concrete object.

Currently, the following feSolutions, fe Parameter and affordances are defined:

{{ inputcpp('src/include/ikarus/finiteElements/feRequirements.hh',False,14,50) }}

## Usage

Inside the finite element the information can than be conveniently extracted:
```cpp linenums="1"
const auto& d      = req.getSolution(FESolutions::displacement);
const auto& lambda = req.getParameter(FEParameter::loadfactor);
if(req.hasAffordance(stiffness))
  ...
```
 and with this you can develop your local finite element.

!!! hint "Affordance"
        It is good style to indicate that you can not fulfill an affordance by throwing an appropriate exception!


# FE result requirements
Additionally, to the upper finite element requirements there are result requirements. 
They have the same methods as finite element requirements but add additional ones.
They are used for the `calculateAt` method of [finite elements](finiteElements.md).
They are a way to communicate the requested results to the finite elements.

## Construction
Similar to above the construction is as follow:
```cpp linenums="1"
ResultRequirements resultRequirements = Ikarus::ResultRequirementsBuilder()
        .insertGlobalSolution(FESolutions::displacement, d)
        .insertParameter(FEParameter::loadfactor, lambda)
        .addResultRequest(ResultType::cauchyStress,,ResultType::director).build();
```

The current supported results are

{{ inputcpp('src/include/ikarus/finiteElements/feRequirements.hh',False,51,61) }}


## Usage
To extract the needed information result requirements have the same interface as finite element requirements.
But they allow to query information which hresults should be calculated.

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