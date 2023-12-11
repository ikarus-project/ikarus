# Master (becomes Release v0.4)

<!--
SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: LGPL-3.0-or-later
-->

- Add CODE_OF_CONDUCT.md file
- Refactored Cmake and directory structure  ([#146](https://github.com/ikarus-project/ikarus/pull/146))
- Added comment section to blog post ([511d83](https://github.com/ikarus-project/ikarus/commit/511d83f9e7c474c9b320db5bc9367114ebe2825d))
- Updated license information ([#138](https://github.com/ikarus-project/ikarus/pull/138))
- Added detailed documentation for ikarus-examples ([#140](https://github.com/ikarus-project/ikarus/pull/140))
- Added computation of Cauchy stress in linear elasticity ([#137](https://github.com/ikarus-project/ikarus/pull/137))
- Added greeting and init function for a reasonable default, e.g., log also to
file ([#147](https://github.com/ikarus-project/ikarus/pull/147))
- Added an interface for the material library ([#154](https://github.com/ikarus-project/ikarus/pull/154))
- Added a wrapper for flat and blocked basis ([#157](https://github.com/ikarus-project/ikarus/pull/157))
- Added basic python bindings and infrastructure. Ikarus can now be downloaded
via `pip install pyikarus` ([#157](https://github.com/ikarus-project/ikarus/pull/152))
- Added result evaluators to post-process desired results ([#165](https://github.com/ikarus-project/ikarus/pull/165))
- Added explicit calculation of scalars, vectors, and matrices for `NonLinearElastic`
([#160](https://github.com/ikarus-project/ikarus/pull/160))
- Renamed Ikarus::linearAlgebraFunctions to Ikarus::functions ([#171](https://github.com/ikarus-project/ikarus/pull/171))
- Added Kirchhoff-Love Shell based on automatic differentiation ([#177](https://github.com/ikarus-project/ikarus/pull/177))
- Added getRawMatrix/getRawVector functionality to the assemblers ([#179](https://github.com/ikarus-project/ikarus/pull/179))
- Add Clang 16 support ([#186](https://github.com/ikarus-project/ikarus/pull/176))

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
- Removed examples from the Ikarus repository; the examples folder is now a sandbox where users can simply compile their own code [#99](https://github.com/ikarus-project/ikarus/pull/99)
- Added class `DirichletValues` to handle homogeneous and inhomogeneous dirichlet values [#104](https://github.com/ikarus-project/ikarus/pull/104)
- Added documentation for ikarus-examples [#106](https://github.com/ikarus-project/ikarus/pull/106)
- Refactored dune <--> Eigen transformations, thereby having smoother transformations between FieldVector and EigenVector/Matrix [#111](https://github.com/ikarus-project/ikarus/pull/111)
- Docker images with Ikarus installed are created from `main`
- Added license statement to each file [#114](https://github.com/ikarus-project/ikarus/pull/114)
- Moved `localfefunctions` to a separate Dune module [#117](https://github.com/ikarus-project/ikarus/pull/117)
- Refactored the documentation and added installation instructions [#125](https://github.com/ikarus-project/ikarus/pull/125)
- Added a workflow to create a release on [GitHub](https://github.com/ikarus-project/ikarus/releases) and push it to [DaRUS](https://darus.uni-stuttgart.de/dataset.xhtml?persistentId=doi%3A10.18419%2Fdarus-3303&version=DRAFT)
- Added the first blog post of Ikarus [#128](https://github.com/ikarus-project/ikarus/pull/128)

## Release v0.2 (Apollodorus)

- Moved from Catch2 testing environment to `Dune::TestSuite`
- Generalized local function test and added `linearStrains` expression
