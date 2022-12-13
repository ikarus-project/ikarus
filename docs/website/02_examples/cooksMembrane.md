<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Cook's membrane
The Cook's membrane problem adapted from the paper[@cook_improved_1974] is implemented in the example
`iks008_cooksMembrane.cpp`. This problem can be solved not only with
structured meshes provided, but also with unstructured and triangular meshes. The input parameters like material and grid
parameters are read from the file `iks008_cooksMembrane.parset`. The problem can be solved also with the standard planar solid element,
or with enhanced assumed strain elements. For more details on the element technologies, refer the
[documentation](../01_framework/finiteElements.md). `iks008` solves the problem with a set of existing finite elements and compares the
convergence rates. 