# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import debug_info

debug_info.setDebugFlags()

import ikarus as iks
import dune.grid
import dune.functions

import math


def makeGrid():
    lowerLeft = []
    upperRight = []
    elements = []
    for _ in range(2):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    return dune.grid.structuredGrid(lowerLeft, upperRight, elements)


def testDirichletValues():
    grid = makeGrid()
    basis = iks.basis(grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2))
    basis = basis.flat()

    dirichletValues = iks.dirichletValues(basis)

    assert basis.size() == dirichletValues.size

    def fixOneIndex(vec, globalIndex):
        if globalIndex == 0:
            vec[globalIndex] = True

    dirichletValues.fixBoundaryDOFs(fixOneIndex)

    # This is equivalent to values[0] == True, but this syntax is discouraged by PEP
    assert dirichletValues.container[0]

    # Note that the result of localView.index(localIndex) is a multiIndex even for a flat basis, the localIndex appears to be a int
    def fixAnotherIndexWithLocalView(vec, localIndex, localView):
        if localView.index(localIndex) == [4]:
            vec[localView.index(localIndex)] = True

        assert isinstance(localIndex, int)
        assert isinstance(localView.index(localIndex), list)
        assert isinstance(localView.index(localIndex)[0], int)

    dirichletValues.fixBoundaryDOFs(fixAnotherIndexWithLocalView)
    assert sum(dirichletValues.container) == 2

    def fixTopSide(vec, localIndex, localView, intersection):
        if intersection.geometry.center[0] == 1.0:
            vec[localView.index(localIndex)] = True

        assert isinstance(localIndex, int)
        assert isinstance(localView.index(localIndex), list)
        assert isinstance(localView.index(localIndex)[0], int)

    dirichletValues.fixBoundaryDOFs(fixTopSide)

    # This assmues a structured grid
    indicesPerDirection: int = (math.sqrt(grid.size(0)) + 1) * 2
    assert sum(dirichletValues.container) == 2 + indicesPerDirection

    ssb0 = dune.functions.subspaceBasis(basis, 0)
    dirichletValues.fixBoundaryDOFsOfSubSpaceBasis(fixTopSide, ssb0)

if __name__ == "__main__":
    testDirichletValues()
