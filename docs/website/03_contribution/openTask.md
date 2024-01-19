# Open tasks

Thank you for your interest in contributing to this code base.
The following task are open for your contributions.

## Local functions

* Implementing a unit normal field function[^KLNote] and its derivatives w.r.t. its coefficients \( \boldsymbol{x}_i \)

    \[
    \boldsymbol{n} = \frac{\boldsymbol{a}_1 \times \boldsymbol{a}_2}{||\boldsymbol{a}_1 \times \boldsymbol{a}_2||}, \quad \text{with }
    \boldsymbol{a}_{\alpha} = \sum_{i=1}^n N^i_{,\alpha}(\boldsymbol{\xi}) \boldsymbol{x}_i
    \]

    To implement this, see [link](../01_framework/localFunctions.md#how-to-implement-your-own-local-functions).

* Support second derivatives
* Add \( \operatorname{div} \) and \( \operatorname{curl} \) wrapper

## Control routines

* Dynamics (Explicit/ implicit time stepping)

## Control routines - addons

* Extended systems
* Nonlinear dependence of $F_{ext}$ on $\mathbf{D}$ and $\lambda$ for path-following techniques, see [control routines](../01_framework/controlRoutines.md#path-following-techniques).

## Finite element helper

* Implement a default mass matrix

## Finite elements

* Nonlinear Reissner-Mindlin shell [@muller2022consistent]
* 3D-Beam
* Standard beam and plate formulations
* Add calculation of different stress measures (push-forward or pull-back of PK2 stresses or similar)

## Further addons

* [Muesli](https://materials.imdea.org/muesli/)

[^KLNote]: This is usually needed for a Kirchhoff-Love shell implementation, see [@kiendlKLshell].

!!! note  "Code style"
    For details on the code style, refer [link](codeStyle.md).

\bibliography
