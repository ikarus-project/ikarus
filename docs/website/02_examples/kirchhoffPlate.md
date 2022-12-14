<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Plate subjected to a surface load
Kirchhoff type plate element is implemented in `iks004_kirchhoffPlate.cpp` using the automatic differentiation
technique as commented before. The basis used for discretization is a NURBS basis from the `dune-iga` module.
The problem is solved and convergence plots are created by comparing the solutions to available analytical solutions for
simply supported and clamped boundaries.