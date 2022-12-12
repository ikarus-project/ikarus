<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Compute the value of $\pi$
The examples `iks002_compute-pi.cpp` and `iks003_compute-pi.cpp`shows the calculation of $\pi$ by computing the area
and circumference of a unit circle. These examples help to understand the `Grid` module from Dune and the refinement techniques it
brings. The example `iks002` shows that a global refinement doesn't refine the number of grid entities on the boundary
of the circle, which leads to a poor approximation of $\pi$ while comparing with the circumference of the circle.
On the other hand, `iks003` shows how elements on the boundaries can be marked and refined, thereby resulting in an
accurate approximation of $\pi$.