<!--
SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
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
- Refactor of FE requirements and non-linear operator construction ([#289](https://github.com/ikarus-project/ikarus/pull/289))
    - The FE requirement type is no longer standalone but is rather dictated by the element.
      Therefore, the elements no longer contain the template argument `FERequirements`.
      It can be obtained statically by `FEType::Requirement` or by getting an object via `fe.createRequirement()`.
    - The exported type by the finite elements is `Requirement` and no longer `FERequirementType`.
    - FE requirements no longer contain affordance information; they now have to be passed separately while calling methods of the assembler.
    - The FE requirements only contain a single solution vector and a single parameter.
    - The getter functions are renamed from `getGlobalSolution` and `getParameter` to `globalSolution` and `parameter`, respectively.
    - The `FErequirements` now have a method `populated` to indicate if the quantities needed by the FE are inserted.
    - The assemblers can now be bound to specific FE requirements, affordances, and an DBCOption.
    - The assemblers now implement a single function `matrix` to assemble matrix quantities.
      Similarly, a single function `vector` is implemented to assemble the corresponding vector quantities.
      Here also, affordance and the `DBCOption` have to be passed.
      The `DBCOption` decides whether a raw, a reduced, or a full matrix is to be obtained.
    - The finite element functions `calculateMatrix`, `calculateVector`, and `calculateScalar` now directly accept the affordances.
    - The affordances are now all singular, `ScalarAffordances` --> `ScalarAffordance`
    - The affordances now know about each other.
      Therefore, calling the function `vectorAffordance(MatrixAffordance::stiffness)` returns `VectorAffordance::forces`.
    - Nonlinear operators can now conveniently be constructed by `NonLinearOperatorFactory`. A bound assembler can be passed to the `op` function,
      but all options can also be passed directly.
    - The linear solver is now fully copyable and moveable.
    - The non-linear solver creation can now be simply done by using the new class `NonlinearSolverFactory`, which receives the settings of
      the non-linear solver and
      adds the function `create`, which accepts a *bound* assembler, where a non-linear operator can be constructed on the fly. For this to
      work, all non-linear solvers have to implement a method `createNonlinearSolver`, which accepts the settings of the class and a nonlinear
      operator. See, e.g., `NewtonRaphsonWithSubsidiaryFunction`.
    - Python:
        - Bindings for Enums can now be done conveniently with the `ENUM_BINDINGS` macro.
        - The finite element functions `calculateMatrix`, `calculateVector`, and `calculateScalar` now directly accept the affordances.
        - The assembler bindings now also accept affordances and `DBCOption`, and they are also renamed  to simply `matrix`, `vector`
          and `scalar`.
        - The assemblers also export the binding functions to bind the assemblers.
- Added a new class `AssemblerManipulator` that wraps an existing assembler
  and helps to manipulate the assembled quantities.
  ([#304](https://github.com/ikarus-project/ikarus/pull/304))
    - This can be used, for instance, to apply concentrated forces or to add spring stiffness in a particular direction.
    - Furthermore, a helper function to get the global index of a Lagrange node at the given global position is added.
- Rework the Python Interface for `DirichletValues` plus add support to easily fix boundary DOFs of `Subspacebasis` in C++ and Python ([#305](https://github.com/ikarus-project/ikarus/pull/305))
- Add a truss element ([#302](https://github.com/ikarus-project/ikarus/pull/302))
- Add an About Ikarus page in the documentation ([#291](https://github.com/ikarus-project/ikarus/pull/291))
- Add new class `Vtk::Writer`, which implements some convenience methods over the existing `dune-vtk` module ([#309](https://github.com/ikarus-project/ikarus/pull/309))
- Add `VanishingStrain` material (useful for example for plane strain case), also refactor the constructor of `LinearElastic` to take any linear material law ([#317](https://github.com/ikarus-project/ikarus/pull/317))
- Add new result evaluators and support for full stress evaluation for the case of vanishing strain and stress in `LinearElastic`
  and `NonlinearElastic` ([#343](https://github.com/ikarus-project/ikarus/pull/343)). Further interface changes:
    - vanishing stress and strain now have a public `underlying()` functionality.
    - the signature of the `operator()` of the result evaluator changed to also include the finite element and the evaluation position.
    - `UserFunction` are now passed directly to the `ResultFunction` and not as a template argument.
    - `Vtk::Writer` now also accepts a `UserFunction` in the `addResult()` function.
- Add hyperelastic and an experimental AutoDiff-based material models ([#333](https://github.com/ikarus-project/ikarus/pull/333)
  and [#342](https://github.com/ikarus-project/ikarus/pull/342))
    - The `Hyperelastic` class takes in its deviatoric and volumetric parts separately as arguments.
        - This class serves as a general interface for hyperelastic material models.
        - It also performs the necessary transformations from principal to Cartesian coordinate systems,
          thereby enabling the deviatoric and volumetric parts to only implement the energy and its derivatives.
    - The `Deviatoric` and `Volumetric` classes in turn accept a deviatoric and a volumetric function, respectively.
        - They serve as an interface that processes the derivatives implemented by the respective functions and
          forwards it to the `Hyperelastic` class.
        - The different deviatoric functions implemented are Arruda-Boyce, Blatz-Ko, Gent, Invariant-based,
          and Ogden models.
        - The invariant-based model is a general model and can be used to create other material models,
          for example, Mooney-Rivlin and Yeoh models.
        - Eleven different volumetric functions are included.
    - An AutoDiff-based material model is included mainly to test these hyperelastic material models. It can be found in the `Experimental` namespace.
    - All materials are now in a separate namespace, `Ikarus::Materials`.
- The solution vector and parameter are now stored internally within the FE requirements instead of being passed by reference.
  This change ensures better encapsulation and simplifies the usage of finite elements. ([#356](https://github.com/ikarus-project/ikarus/pull/356))
- Add wrapper classes for small strain and finite strain materials from the Muesli material library
  ([#337](https://github.com/ikarus-project/ikarus/pull/337))
- `Observers` and `Observables` are replaced with `Broadcasters` and `Listeners`. Existing loggers work almost the same.
 A noteworthy difference is that, the Broadcaster (e.g. a nonlinear solver) has to be registered to a Listener (e.g. a logger) with `logger.subscribeTo(solver)` ([#349](https://github.com/ikarus-project/ikarus/pull/349))
- Refactor EAS to handle NonLinearElastic ([#325](https://github.com/ikarus-project/ikarus/pull/325))
    - This also includes a new interface method in the `Mixin` called `updateState` with the respective `impl()` function.
    - For EAS, `updateStateImpl()` is used to update the internal variable `alpha_` in a nonlinear analysis.
    - Missing functions like `getStress` and `materialTangentFunction` are added to `LinearElastic` and `NonLinearElastic`, respectively.
    - `EnhancedAssumedStrains` class now takes a struct denoting the enhanced strain type as a template argument, which implements the
      enhanced strain value, its first and second derivatives w.r.t the displacements and the internal variable `alpha_`.
- Add Q1Hn and Q1HTn elements ([#363](https://github.com/ikarus-project/ikarus/pull/363))
- Add wrappers for solving the generalized eigenvalue problem with sparse and dense matrices with Eigen and Spectra ([#368](https://github.com/ikarus-project/ikarus/pull/368)).
- Add functionality to invert hyperelastic material laws. For Neo-Hooke and SVK law, analytical solutions are available.
  For the general hyperelastic framework, a numerical approach is used ([#369](https://github.com/ikarus-project/ikarus/pull/369)).
- Add `AssumedStress` class to extend linear and nonlinear solid elements. The implementation is similar to EAS elements, hence the
  internal variables here and in EAS were renamed to `internalVariable` ([#371](https://github.com/ikarus-project/ikarus/pull/371)).
- Add push-forward and pull-back operations for stress measures ([#372](https://github.com/ikarus-project/ikarus/pull/372)).

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
- Add Clang 16 support ([#186](https://github.com/ikarus-project/ikarus/pull/186))
- Added a default adaptive step sizing possibility and refactored loggers ([#193](https://github.com/ikarus-project/ikarus/pull/193))
- Added a wrapper to fix Dirichlet BCs for Lagrange Nodes ([#222](https://github.com/ikarus-project/ikarus/pull/222))
- Refactored `getDisplacementFunction` in finite elements ([#223](https://github.com/ikarus-project/ikarus/pull/223))
- Refactored load functions in finite elements ([#221](https://github.com/ikarus-project/ikarus/pull/221))
- Added doxygen class documentation ([#220](https://github.com/ikarus-project/ikarus/pull/220))

- Improve material library and Python bindings ([#186](https://github.com/ikarus-project/ikarus/pull/186)), default e.g.
  `StVenantKirchhoff` is
  not a template anymore use `StVenantKirchhoffT<...>` scalar-types that are not `double`

- `Ikarus::ILinearSolver` is no template anymore and renamed and uses double, use
  `Ikarus::LinearSolverT<...>` instead ([#186](https://github.com/ikarus-project/ikarus/pull/186)),
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
