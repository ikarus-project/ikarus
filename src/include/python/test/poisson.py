import numpy as np
from scipy.sparse import lil_matrix

import dune.geometry
import dune.grid
import dune.functions

# f(x) = 1
f = lambda x:  1

# g(x) = 0
g = lambda x:  0

# boundary of the unit square
def isDirichlet(x):
    return x[0] <= 0.0 or x[0] >= 1.0 or x[1] >= 1.0 or x[1] <= 0.0


# TODO: This assembler loop is very inefficient in terms of run time and should be improved using Python vectorization.
# See discussion at https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/295 for hints and code pointers.
def localAssembler(element,localBasis):

    n = len(localBasis)

    localA = np.zeros((n,n))
    localB = np.zeros(n)

    # choose a high enough quadrature order
    quadOrder = 3

    # create a quadrature rule and integrate
    quadRule = dune.geometry.quadratureRule(element.type, quadOrder)
    for pt in quadRule:

        pos = pt.position

        # evaluate the local basis functions (this is an array!)
        phi = localBasis.evaluateFunction(pos)

        # evaluate the local basis reference Jacobians (array of arrays)
        phiRefGradient = localBasis.evaluateJacobian(pos)

        # get the transformation matrix
        jacobianInverseTransposed = element.geometry.jacobianInverseTransposed(pos)

        # |det(jacobianInverseTransposed)|
        integrationElement = element.geometry.integrationElement(pos)

        # transform the reference Jacobians to global geometry
        phiGradient = [ np.dot(jacobianInverseTransposed, np.array(g)[0]) for g in phiRefGradient ]

        posGlobal = element.geometry.toGlobal(pos)

        for i in range( n ):
            for j in range( n ):
                # compute grad(phi_i) * grad(phi_j)
                localA[i,j] += pt.weight * integrationElement * np.dot( phiGradient[i], phiGradient[j] )

            localB[i] += pt.weight * integrationElement * phi[i] * f(posGlobal)

    return localA, localB



def globalAssembler(basis):

    N = len(basis)

    A = lil_matrix( (N,N) )
    b = np.zeros( N )

    # mark all Dirichlet DOFs
    dichichletDOFs = np.zeros( N, dtype=bool )
    basis.interpolate(dichichletDOFs,isDirichlet)

    # interpolate the boundary values
    gCoeffs = np.zeros( N )
    basis.interpolate(gCoeffs,g)

    # extract grid and localView
    localView = basis.localView()
    grid = basis.gridView

    for element in grid.elements:

        # assign the localView to the current element
        localView.bind(element)

        # set of all shape functions with support in this element
        localBasis = localView.tree().finiteElement.localBasis

        localN = len(localBasis)

        localA, localb = localAssembler(element,localBasis)

        # copy the local entries into the global matrix using the
        # index mapping given by the localView
        for i in range(localN):

            gi = localView.index(i)[0]

            if dichichletDOFs[gi]:
                A[gi,gi] = 1.0
                b[gi] = gCoeffs[gi]
            else:
                b[gi] += localb[i]
                for j in range(localN):
                    gj = localView.index(j)[0]
                    A[gi, gj] += localA[i, j]

    return A,b


############################### START ########################################

# number of grid elements (in one direction)
gridSize = 4

# create a grid of the unit square
grid = dune.grid.structuredGrid([0,0],[1,1],[gridSize,gridSize])

# create a nodal Lagrange FE basis of order 1
basis = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=1))

# compute A and b
A,b = globalAssembler(basis)

# convert A to a dense matrix for testing
A_dense = A.todense()

# solve!
x = np.linalg.solve(A_dense, b)

# test the result
xTest = [ 0., 0., 0., 0., 0., 0., 0.04821429, 0.06026786, 0.04821429, 0.,
        0., 0.06026786, 0.07767857, 0.06026786, 0., 0., 0.04821429,
        0.06026786, 0.04821429, 0., 0., 0., 0., 0., 0. ]

assert (np.linalg.norm(x-xTest) < 1e-7)
