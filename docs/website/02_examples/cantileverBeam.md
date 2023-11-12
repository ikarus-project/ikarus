# Cantilever beam with point load

## Description

The example `iks002_cantileverBeamOneDGrid.cpp` demonstrates a simple implementation of the standard one-dimensional
Timoshenko beam element, which is clamped on the left-hand side. A point load is applied on the right-hand side of the beam. It uses
`Dune::OneDGrid` to generate the required
grid. The implementation shown here assembles the stiffness matrices explicitly. Advanced
implementations of matrix assembly and other features of Ikarus are showcased in the other examples.

## Code highlights

The `#!cpp main()` function calls the following function to run the example.

```cpp
void exampleTimoshenkoBeam(const int polynomialOrderW, // (1)!
const int polynomialOrderPhi, // (2)!
const int numElements // (3)!
);
```

1. polynomial order for the approximation of the displacement `w`
2. polynomial order for the approximation of the rotation `phi`
3. number of elements

This function has pre-defined geometric and material parameters that are mentioned below:

```cpp
const double b  = 1; // (1)!
const double L  = 10; // (2)!
const double E  = 1000; // (3)!
const double G  = E / 2;  // Poisson's ratio = 0 // (4)!
const double t  = 1e-3; // (5)!
```

1. breadth of the beam (rectangular cross-section assumed)
2. length of the beam
3. Young's modulus
4. Shear modulus (Poisson's ratio is assumed to be zero)
5. thickness of the beam (rectangular cross-section assumed)

The material matrix is then defined as

```cpp
Eigen::Matrix2d C;
C << EI, 0, 0, GA;
```

The one-dimension grid is created by specifying the starting (`0`) and end (`L`) points of the beam and using `numElements`.
The grid module `Dune::OneDGrid`[@sander2020dune] is used here. The *composite basis* feature in Dune allows for
different bases for the degrees of freedom `w` and `phi`. Empty dense global and local stiffness matrices are then defined along with
a global external load vector. The quadrature rule for integration of the stiffness matrix and the external load vector is set up as shown below:

```cpp
const int maxOrderIntegration = std::max(2 * (polynomialOrderW - 1), 2 * polynomialOrderPhi);
const auto &rule
    = Dune::QuadratureRules<double, 1>::rule(ele.type(), maxOrderIntegration, Dune::QuadratureType::GaussLegendre);
```

The local element stiffness matrix is then obtained by the following function:

```cpp
void TimoshenkoBeamStiffness(auto &KLocal, auto &localView, auto &gridElement, auto &quadratureRule, const Eigen::Matrix2d &C);
```

Access to the information corresponding to the two independent degrees of freedom (`w` and `phi`) can be obtained by using the child
elements as depicted below:

```cpp
auto wFE   = localView.tree().child(_0);
auto phiFE = localView.tree().child(_1);
Dune::CachedLocalBasis basisW(wFE.finiteElement().localBasis()); // (1)!
Dune::CachedLocalBasis basisPhi(phiFE.finiteElement().localBasis()); // (2)!
```

1. a local basis created for `w`
2. a local basis created for `phi`

The determinant of Jacobian is obtained by using `#!cpp auto detJ = gridElement.geometry().volume();`. The
$\mathbf{B}$-operator and the local stiffness matrix are then obtained by

```cpp
Eigen::Matrix2Xd B;
for (auto &gp : quadratureRule) { // (1)!
  basisW.evaluateJacobian(gp.position(), dNwDxi); // (2)!
  basisPhi.evaluateFunction(gp.position(), Nphi); // (3)!
  basisPhi.evaluateJacobian(gp.position(), dNphiDxi); // (4)!

  B.setZero(Eigen::NoChange, numDofsPerEle);

  for (unsigned int i = 0; i < wFE.size(); ++i)
    B(1, wFE.localIndex(i)) = dNwDxi[i] / detJ; // (5)!

  for (unsigned int i = 0; i < phiFE.size(); ++i)
    B.col(phiFE.localIndex(i)) << dNphiDxi[i] / detJ, Nphi[i]; // (6)!

  KLocal += B.transpose() * C * B * detJ * gp.weight(); // (7)!
}
```

1. Looping over the integration point
2. derivative of the ansatz function of `w` with respect to the local parametric space
3. ansatz function of `phi`
4. derivative of the basis of `phi` with respect to the local parametric space
5. Filling up the $\mathbf{B}$-operator for the columns corresponding to 'w'
6. Filling up the $\mathbf{B}$-operator for the columns corresponding to 'phi'
7. Integrating to arrive at the local stiffness matrix

Assembly of the local element stiffness matrices is done to obtain the global element stiffness matrix, as shown below:

```cpp
for (auto i = 0U; i < localView.size(); ++i)
  for (auto j = 0U; j < localView.size(); ++j)
    KGlobal(localView.index(i)[0], localView.index(j)[0]) += KLocal(i, j);
```

The point load on the right end of the beam is applied by setting the corresponding entry in the global external load
vector to the prescribed value `F` using the command `#!cpp FExtGlobal(getGlobalDofId(TimoshenkoBeam::w, basis, L)) = F;`.
The function

```cpp
unsigned int getGlobalDofId(TimoshenkoBeam requestedQuantity, const auto &basis, const double position);
```

with

```cpp
enum class TimoshenkoBeam { w, phi };
```

is used to get the degree of freedom of the `requestedQuantity` (`w` or `phi`) for the defined `basis` at position `L`.

The left end of the beam is clamped by using the following code:

```cpp
std::vector<unsigned int> fixedDofs{getGlobalDofId(TimoshenkoBeam::w, basis, 0.0),
                                    getGlobalDofId(TimoshenkoBeam::phi, basis, 0.0)};
for (auto dof : fixedDofs) {
  KGlobal.col(dof).setZero();
  KGlobal.row(dof).setZero();
  KGlobal(dof, dof) = 1.0;
}
```

Finally, the system of equations is solved by using the solver type `#!cpp Ikarus::SolverTypeTag::d_LDLT`. For more
details on the solver types, refer to the [documentation](../01_framework/solvers.md).
For post-processing, the deformed configuration of the beam can be plotted using the following function, shown here in the example:

```cpp
void plotDeformedTimoschenkoBeam(auto &gridView, auto &basis, auto &d_glob, double EI, double GA, double L, double F);
```

This function uses the plotting features of [Matplot++](https://github.com/alandefreitas/matplotplusplus),
which has a similar syntax to [Matplotlib](https://matplotlib.org/).

## Takeaways

- `#!cpp Dune::OneDGrid` can be used to create one-dimensional grids.
- Grids and quadrature rules from Dune can be directly incorporated into the Ikarus framework.
- A composite basis can be used to have different ansatz functions for different degrees of freedom.
- Solvers from the Eigen library can be used to solve the linear system of equations.
- `localBasis` functions can be used to evaluate the ansatz functions and its derivatives.
