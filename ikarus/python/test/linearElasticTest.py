# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import setpath

setpath.set_path()
import ikarus as iks
import ikarus.finite_elements
import ikarus.utils
import ikarus.assembler
import ikarus.dirichletValues
import numpy as np
import scipy as sp

import dune.grid
import dune.functions
from dune.vtk import vtkWriter

if __name__ == "__main__":
    lowerLeft = []
    upperRight = []
    elements = []
    for i in range(2):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    req = ikarus.FErequirements()
    req.addAffordance(iks.ScalarAffordances.mechanicalPotentialEnergy)

    grid = dune.grid.structuredGrid(lowerLeft, upperRight, elements)
    grid.hierarchicalGrid.globalRefine(6)
    basisLagrange12 = dune.functions.defaultGlobalBasis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    )
    dirichletValuesT = iks.dirichletValues(basisLagrange12)

    basisLagrange1 = iks.basis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    )
    flatBasis = basisLagrange1.flat()
    d = np.zeros(len(flatBasis))
    d[0] = 0.0

    lambdaLoad = iks.ValueWrapper(3.0)
    req.insertParameter(iks.FEParameter.loadfactor, lambdaLoad)

    assert req.getParameter(iks.FEParameter.loadfactor) == lambdaLoad
    req.insertGlobalSolution(iks.FESolutions.displacement, d)

    d2 = req.getGlobalSolution(iks.FESolutions.displacement)

    assert ("{}".format(hex(d2.__array_interface__["data"][0]))) == (
        "{}".format(hex(d.__array_interface__["data"][0]))
    )
    assert len(d2) == len(d)
    assert (d2 == d).all()
    fes = []

    def volumeLoad(x, lambdaVal):
        return np.array([lambdaVal * x[0] * 2, 2 * lambdaVal * x[1] * 0])

    def neumannLoad(x, lambdaVal):
        return np.array([lambdaVal * 0, lambdaVal])

    neumannVertices = np.zeros(grid.size(2) * 2, dtype=bool)
    assert len(neumannVertices) == len(flatBasis)

    flatBasis.interpolate(neumannVertices, lambda x: True if x[1] > 0.9 else False)

    boundaryPatch = iks.utils.boundaryPatch(grid, neumannVertices)

    # the following should throw
    try:
        for e in grid.elements:
            iks.finite_elements.linearElasticElement(
                basisLagrange1, e, 1000, 0.2, volumeLoad, boundaryPatch
            )
        assert False
    except TypeError:
        pass

    for e in grid.elements:
        fes.append(
            iks.finite_elements.linearElasticElement(
                basisLagrange1, e, 1000, 0.2, volumeLoad, boundaryPatch, neumannLoad
            )
        )

    forces = np.zeros(8)
    stiffness = np.zeros((8, 8))
    fes[0].calculateVector(req, forces)
    fes[0].calculateMatrix(req, stiffness)
    fes[0].localView()

    dirichletValues = iks.dirichletValues(flatBasis)

    def fixFirstIndex(vec, globalIndex):
        vec[0] = True

    def fixAnotherVertex(vec, localIndex, localView):
        localView.index(localIndex)
        vec[1] = True

    def fixLeftHandEdge(vec, localIndex, localView, intersection):
        if intersection.geometry.center[1] < -0.9:
            vec[localView.index(localIndex)] = True

    dirichletValues.fixBoundaryDOFs(fixFirstIndex)
    dirichletValues.fixBoundaryDOFsUsingLocalView(fixAnotherVertex)
    dirichletValues.fixBoundaryDOFsUsingLocalViewAndIntersection(fixLeftHandEdge)

    assembler = iks.assembler.sparseFlatAssembler(fes, dirichletValues)
    assemblerDense = iks.assembler.denseFlatAssembler(fes, dirichletValues)

    Msparse = assembler.getMatrix(req)
    forces = assembler.getVector(req)

    x = sp.sparse.linalg.spsolve(Msparse, -forces)
    fx = flatBasis.asFunction(x)
    grid.plot()

    writer = vtkWriter(grid, "nameTest", pointData={("displacement", (0, 1)): fx})
