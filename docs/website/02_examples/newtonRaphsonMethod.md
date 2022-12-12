<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Newton-Raphson method
`iks006_newtonRaphson.cpp` shows a basic example of the Newton-Raphson method to solve a non-linear set of equations.
A function which shows the algorithm explicitly is provided and another function which is implemented in Ikarus is
demonstrated. The function which depicts the Ikarus implementation uses a
[non-linear operator](https://ikarus-project.github.io/01_theory/nonlinearOperator/) to
perform the Newton-Raphson iterations. A logger can also be subscribed to in order to observe the residual norms,
for instance.