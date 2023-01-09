import dune.grid
import dune.functions


# Test all currently supported global bases in dimension `dimension`.
# Only the size of the global bases are tested.
def test(dimension):
    lowerLeft = []
    upperRight = []
    elements = []
    for i in range(dimension):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    grid = dune.grid.structuredGrid(lowerLeft,upperRight,elements)
    nDofs1 = grid.size(dimension)
    nDofs2 = grid.size(dimension)
    for i in range(dimension):
        nDofs2 += grid.size(i)

    basisLagrange1 = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=1))
    assert(len(basisLagrange1) == nDofs1)

    basisLagrange2 = dune.functions.defaultGlobalBasis(grid, dune.functions.Lagrange(order=2))
    assert(len(basisLagrange2) == nDofs2)

    basisPower = dune.functions.defaultGlobalBasis(grid, dune.functions.Power(dune.functions.Lagrange(order=1),exponent=dimension))
    assert(len(basisPower) == nDofs1*dimension)

    basisComposite = dune.functions.defaultGlobalBasis(grid, dune.functions.Composite(dune.functions.Power(dune.functions.Lagrange(order=2),exponent=dimension),
                                                                                      dune.functions.Lagrange(order=1)))
    assert(len(basisComposite) == nDofs2*dimension + nDofs1)
