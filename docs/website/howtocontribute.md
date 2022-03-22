# How to contribute
Thanks for your interest in contributing to this code base.
If your are interested the following task are vacant


### Local functions

Implementing a unit normal field function[^KLNote] and its derivatives w.r.t. its coefficients $\boldsymbol{x}_i$

$$ 
\boldsymbol{n} = \frac{\boldsymbol{a}_1 \times \boldsymbol{a}_2}{||\boldsymbol{a}_1 \times \boldsymbol{a}_2||}, \quad \text{with } \boldsymbol{a}_{\alpha} = \sum_{i=1}^n N^i_{,\alpha}(\boldsymbol{\xi}) \boldsymbol{x}_i
$$ 


To implement these see [link](01_theory/LocalFunctions.md#how-to-implement-your-own-local-functions)


[^KLNote]: This is usually needed for a Kirchhoff-Love shell implementation, see [@kiendlKLshell].

\bibliography 