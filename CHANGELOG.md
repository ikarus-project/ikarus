<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: LGPL-2.1-or-later
-->

# Master (becomes Release v0.3)

- Added linear solver capability to accept matrix-valued rhs #76
- Added `EnhancedAssumedStrains` that can be used to decorate a linear elasticity element with several EAS methods #74
- Added codespell workflow #70 (CI checks now for grammar and typos in comments and variable names)
- Added Path following technique to using, e.g. `Arc length method` the scalar subsidiary equation can be implemented independently #80
- Separated examples from Ikarus repository, now the examples folder is a sandbox to simple compile some user code  #99 +
- Refactored dune <--> Eigen transformations, now having a FieldVector or EigenVector/Matrix, the transformations should feel smoother #111
- Docker container for main are automatically generated to test the examples, #52c4d859e721479a218575b53a7a599df762fdc0
- Added documentation for the examples, #ebfe2f353778f6da9da5aa64be498eead2c3490c
- Added class `DirichletValues` to take care of homogeneous and inhomogeneous dirichlet Values #104
- Added license statement to each file [#114](https://github.com/IkarusRepo/Ikarus/pull/114)

## Release v0.2 (Apollodorus)

- Moved from Catch2 testing environment to `Dune::TestSuite`
- Generalized local function test and added `linearStrains` expression

