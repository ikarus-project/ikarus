<!--
SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: LGPL-3.0-or-later
-->

# Master (becomes Release v0.5)

- Add tensorProduct quadrature rule and [0,1] to [-1,1] transformation ([#238](https://github.com/ikarus-project/ikarus/pull/238))
- Refactored Interface for `calculateAt`-functions and `ResultFunction` by getting rid of `ResultRequirements`
  ([#245](https://github.com/ikarus-project/ikarus/pull/245))
- Refactor `enhancedassumedstrains.hh` ([#247](https://github.com/ikarus-project/ikarus/pull/247))
- add factory functions to `ResultFunction` ([#251](https://github.com/ikarus-project/ikarus/pull/251))
- Refactor `febases` to handle any type of basis ([#256](https://github.com/ikarus-project/ikarus/pull/256))
- Move `autodiffhelper.hh` from `utils` to `autodiff` folder ([#259](https://github.com/ikarus-project/ikarus/pull/259))
- Added a unified, generic interface for `ResultTypes` with the ability to register new
  ones ([#261](https://github.com/ikarus-project/ikarus/pull/261))
- Big refactor of fe requirements and non linear operator construction ([#289](https://github.com/ikarus-project/ikarus/pull/289))
  - The fe requirement type is not any longer standalone but it is dictated by the element. Therefore, the elements are no longer a template of the requirements. It can be obtained statically by `FEType::Requirement` or obtaining an object be `fe.createRequirement()`
  - The exported type by the finite elements is `Requirement` and no longer `FERequirementType`
  - FE requirements no longer contain the affordance information, they now have to be passed separate.
  - The fe requirements only contain now a single solution vector and a single parameter.
  - The getter are renamed from `getGlobalSolution` and `getParameter` to without the `get`.
  - The fe requirements now have a method `populated` to indicated, if the needed quantities are inserted.
  - The assemblers can now be bound to specific fe requirements, affordances and an EnforcingDBCOption.
  - The assemblers now only implement a single function `matrix` and `vector`, where also affordance and the `EnforcingDBCOption` have to be passed. The
    EnforcingDBCOption decides whether a raw, a reduced or a full matrix is obtained. Full means that the Dirichlet boundary conditions are written in the matrix but keeps the size
  - The finite element functions `calculateMatrix`,`calculateVector` and `calculateScalar` now directly accept the affordances
  - The affordances are now all singular `ScalarAffordances` --> `ScalarAffordance`
  - The affordances now also now about each others relation therefore calling the function `vectorAffordance(MatrixAffordance::stiffness)` returns `VectorAffordance::forces`
  - Nonlinear operators can now conveniently be constructed by `NonLinearOperatorFactory`. A bound assembler can be passed to the `op` function but all options can also be passed directly
  - The linear solver are now fully copyable and moveable
  - The non-linear solver creation can now be simply done by using the new class `NonlinearSolverFactory` which receives the settings of the non-linear solver and
    adds the function `.create`, which accepts a *bound* assembler, where a non-linear operator can be constructed on the fly. For this to work all non-linear solvers have to implement a method `createNonlinearSolver`, which accepts the settings of the class and a nonlinear operator. See e.g. `NewtonRaphsonWithSubsidiaryFunction`.

  - Python:
    - Bindings for Enums can now be done conveniently with the `ENUM_BINDINGS` macro
    - The finite element functions `calculateMatrix`,`calculateVector` and `calculateScalar` now directly accept the affordances
    - The assembler bindings now also accept affordances and `EnforcingDBCOption` and they are also renamed to simply `matrix`, `vector`, `scalar`
    - The assemblers also export the binding functions to bind the assemblers



## Release v0.4 (Ganymede)

- Add CODE_OF_CONDUCT.md file
- Refactored Cmake and directory structure ([#146](https://github.com/ikarus-project/ikarus/pull/146))
- Added comment section to blog post ([511d83](https://github.com/ikarus-project/ikarus/commit/511d83f9e7c474c9b320db5bc9367114ebe2825d))
- Updated license information ([#138](https://github.com/ikarus-project/ikarus/pull/138))
- Added detailed documentation for ikarus-examples ([#140](https://github.com/ikarus-project/ikarus/pull/140))
- Added computation of Cauchy stress in linear elasticity ([#137](https://github.com/ikarus-project/ikarus/pull/137))
- Added greeting and init function for a reasonable default, e.g., log also to
  file ([#147](https://github.com/ikarus-project/ikarus/pull/147))
- Added an interface for the material library ([#154](https://github.com/ikarus-project/ikarus/pull/154))
- Added a wrapper for flat and blocked basis ([#157](https://github.com/ikarus-project/ikarus/pull/157))
- Added basic python bindings and infrastructure. Ikarus can now be downloaded
  via `pip install pyikarus` ([#152](https://github.com/ikarus-project/ikarus/pull/152))
- Added result evaluators to post-process desired results ([#165](https://github.com/ikarus-project/ikarus/pull/165))
- Added explicit calculation of scalars, vectors, and matrices for `NonLinearElastic`
  ([#160](https://github.com/ikarus-project/ikarus/pull/160))
- Renamed Ikarus::linearAlgebraFunctions to Ikarus::functions ([#171](https://github.com/ikarus-project/ikarus/pull/171))
- Added Kirchhoff-Love Shell based on automatic differentiation ([#177](https://github.com/ikarus-project/ikarus/pull/177))
- Added getRawMatrix/getRawVector functionality to the assemblers ([#179](https://github.com/ikarus-project/ikarus/pull/179))
- Add Clang 16 support ([#186](https://github.com/ikarus-project/ikarus/pull/176))
- Added a default adaptive step sizing possibility and refactored loggers ([#193](https://github.com/ikarus-project/ikarus/pull/193))
- Added a wrapper to fix Dirichlet BCs for Lagrange Nodes ([#222](https://github.com/ikarus-project/ikarus/pull/222))
- Refactored `getDisplacementFunction` in finite elements ([#223](https://github.com/ikarus-project/ikarus/pull/223))
- Refactored load functions in finite elements ([#221](https://github.com/ikarus-project/ikarus/pull/221))
- Added doxygen class documentation ([#220](https://github.com/ikarus-project/ikarus/pull/220))

- Improve material library and Python bindings ([#186](https://github.com/ikarus-project/ikarus/pull/176)), default e.g.
  `StVenantKirchhoff` is
  not a template anymore use `StVenantKirchhoffT<...>` scalar-types that are not `double`

- `Ikarus::ILinearSolver` is no template anymore and renamed and uses double, use
  `Ikarus::LinearSolverT<...>` instead ([#186](https://github.com/ikarus-project/ikarus/pull/176)),
  if you need a fancy scalar type

## Release v0.3 (Prometheus)

- Added codespell workflow (CI checks now for grammar and typos in comments
  and variable names) [#70](https://github.com/ikarus-project/ikarus/pull/70)
- Added `EnhancedAssumedStrains` to decorate a linear-elastic element with
  various EAS methods [#74](https://github.com/ikarus-project/ikarus/pull/74)
- Added the ability for the linear solver to accept matrix-valued rhs [#76](https://github.com/ikarus-project/ikarus/pull/76)
- Added a path-following technique, such that a scalar subsidiary equation, for example, for `Arc length method`,
  can be implemented independently [#80](https://github.com/ikarus-project/ikarus/pull/80)
- Removed examples from the Ikarus repository; the examples folder is now a sandbox where users can simply compile their own
  code [#99](https://github.com/ikarus-project/ikarus/pull/99)
- Added class `DirichletValues` to handle homogeneous and inhomogeneous dirichlet
  values [#104](https://github.com/ikarus-project/ikarus/pull/104)
- Added documentation for ikarus-examples [#106](https://github.com/ikarus-project/ikarus/pull/106)
- Refactored dune <--> Eigen transformations, thereby having smoother transformations between FieldVector and
  EigenVector/Matrix [#111](https://github.com/ikarus-project/ikarus/pull/111)
- Docker images with Ikarus installed are created from `main`
- Added license statement to each file [#114](https://github.com/ikarus-project/ikarus/pull/114)
- Moved `localfefunctions` to a separate Dune module [#117](https://github.com/ikarus-project/ikarus/pull/117)
- Refactored the documentation and added installation instructions [#125](https://github.com/ikarus-project/ikarus/pull/125)
- Added a workflow to create a release on [GitHub](https://github.com/ikarus-project/ikarus/releases) and push it
  to [DaRUS](https://darus.uni-stuttgart.de/dataset.xhtml?persistentId=doi%3A10.18419%2Fdarus-3303&version=DRAFT)
- Added the first blog post of Ikarus [#128](https://github.com/ikarus-project/ikarus/pull/128)

## Release v0.2 (Apollodorus)

- Moved from Catch2 testing environment to `Dune::TestSuite`
- Generalized local function test and added `linearStrains` expression
