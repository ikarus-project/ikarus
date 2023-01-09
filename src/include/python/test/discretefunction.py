from dune.grid import structuredGrid
import dune.functions
import numpy as np

# create a simple grid
grid = structuredGrid([0,0],[1,1],[10,10])

# simple P1 basis
basis = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=1))

# boring coefficient vector
coeff = np.zeros( len(basis) )

# discrete function over the grid
discreteFunction = basis.asFunction(coeff)

# extract the original grid
grid2 = discreteFunction.grid

# check if there are still 100 squares
assert( grid2.size(0) == 100 )
