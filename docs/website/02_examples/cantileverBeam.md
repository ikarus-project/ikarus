<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Cantilever beam with point load
The example `iks002_cantileverBeamOneDGrid.cpp` shows a simple implementation of a one dimensional Timoshenko beam which is clamped on the left 
hand side. A point load is applied on the right hand side of the beam. It uses `Dune::OneDGrid` to generate the required 
grid. A simple implementation is shown here where the stiffness matrices are assembled explicitly. Advanced 
implementations of matrix assembly and other features of Ikarus is showcased in the other examples.