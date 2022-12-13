<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Non-linear Elasticity for 2D solids
Again, automatic differentiation based implementation is used to perform a non-linear analysis for a 2D block in
`iks006_nonlinear2DSolid.cpp`. Various methods to obtain a 2D grid via Dune is also shown in the commented section in
the beginning. Python is used to provide a Neumann boundary condition providing a demonstration for the usage of a
Python-based code within the Ikarus framework. Load control method is chosen as the desired control routine and
Newton-Raphson (or Trust region methods) are used to solve the non-linear problem itself.