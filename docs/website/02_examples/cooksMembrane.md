<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Cook's membrane
The Cook's membrane problem adapted from the paper[@cook_improved_1974] is implemented in examples
`iks009_cook_membrane.cpp` and `iks010_cook_membrane_convergence.cpp`. This problem can be solved not only with
structured meshes provided, but also with unstructured and triangular meshes. The input parameters like material and grid
parameters are read from the file `cook.parset`. The problem can be solved also with the standard planar solid element,
or with enhanced assumed strain elements. For more details on the element technologies, refer the
[documentation](https://ikarus-project.github.io/01_theory/finiteElements/). `iks009` solves the problem for a chosen
finite element type whereas `iks010` solves the problem with a set of existing finite elements and compares the
convergence rates. 