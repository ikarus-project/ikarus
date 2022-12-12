<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Compression of an incompressible rubber block
`iks004_incompressible_LinearElasticity.cpp` uses a finite element technology with displacement and pressure as
independent degrees of freedom to simulate the compression of an incompressible rubber block. The potential energy for such a system is defined in the
`calculateScalarImpl(const FERequirementType &par, const Eigen::VectorX<ScalarType> &dx)` function in `struct`
named `Solid`. This function uses the principles of automatic differentiation to provide the stiffness matrices and
other necessary quantities to provide a static structural analysis.