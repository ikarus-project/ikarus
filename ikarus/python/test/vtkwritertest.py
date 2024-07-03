# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later


import debug_info

debug_info.unsetDebugFlags()

import math

import ikarus as iks
from ikarus import finite_elements, utils, assembler, io
import numpy as np
import scipy as sp

import dune.grid
import dune.functions

from dune.grid import gridFunction
from dune.vtk import FormatTypes, DataTypes


def runTest():
    lowerLeft = []
    upperRight = []
    elements = []
    for i in range(2):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    grid = dune.grid.structuredGrid(lowerLeft, upperRight, elements)

    basisLagrange1 = iks.basis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    )
    flatBasis = basisLagrange1.flat()
    d = np.zeros(len(flatBasis))
    d[0] = 0.0

    lambdaLoad = iks.ValueWrapper(3.0)

    fes = []

    def vL(x, lambdaVal):
        return np.array([lambdaVal * x[0] * 2, 2 * lambdaVal * x[1] * 0])

    vLoad = iks.finite_elements.volumeLoad2D(vL)

    linElastic = iks.finite_elements.linearElastic(youngs_modulus=1000, nu=0.2)
    for e in grid.elements:
        fes.append(iks.finite_elements.makeFE(basisLagrange1, linElastic, vLoad))
        fes[-1].bind(e)

    req = fes[0].createRequirement()
    req.insertParameter(lambdaLoad)
    req.insertGlobalSolution(d)

    dirichletValues = iks.dirichletValues(flatBasis)

    def fixLeftHandEdge(vec, localIndex, localView, intersection):
        if intersection.geometry.center[1] < -0.99999:
            vec[localView.index(localIndex)] = True

    dirichletValues.fixBoundaryDOFsUsingLocalViewAndIntersection(fixLeftHandEdge)

    sparseAssembler = iks.assembler.sparseFlatAssembler(fes, dirichletValues)
    sparseAssembler.bind(
        req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full
    )

    Msparse = sparseAssembler.matrix()
    forces = sparseAssembler.vector()

    x = sp.sparse.linalg.spsolve(Msparse, -forces)
    req.insertGlobalSolution(x)

    _ = iks.io.vtkWriter(
        sparseAssembler, datatype=DataTypes.Float64, headertype=DataTypes.UInt32
    )


    writer = iks.io.vtkWriter(sparseAssembler, format=FormatTypes.ascii)
    writer.addInterpolation(x, flatBasis, "displacements", 2)


    writer.addAllResultsAsCellData()
    writer.addAllResultsAsPointData()

    writer.write("file")

    assert iks.io.dataCollectors[0] == "lagrange"
    assert iks.io.dataCollectors[1] == "discontinous"

    writer2 = iks.io.vtkWriter(sparseAssembler, dataCollector="lagrange", order=2)

    writer2.addResultAsCellData("linearStress")
    writer2.addResultAsPointData("linearStress")

    writer2.setFormat(FormatTypes.ascii)
    writer2.setDatatype(DataTypes.Float64)
    writer2.setHeadertype(DataTypes.UInt16)

    @gridFunction(grid)
    def g(x):
        return [math.sin(2 * math.pi * x[0] * x[1]), x[0] * x[1]] * 5

    writer2.addPointData(g, name="g", components=(0, 1))
    writer2.addPointData(g, name="g2", components=[0, 1, 2])

    writer2.write("file2")

    # Structured writer
    writer3 = iks.io.vtkWriter(
        sparseAssembler, structured=True, format=FormatTypes.ascii
    )
    writer3.addAllResultsAsCellData()
    fileName = writer3.write("file3")
    assert fileName[-3:] == "vtr"

    # These two constructors should result in a warning
    _ = iks.io.vtkWriter(sparseAssembler, dataCollector="lagrange", structured=True)
    _ = iks.io.vtkWriter(sparseAssembler, dataCollector="unknown")

    
    discontinousVtkWriter = iks.io.vtkWriter(sparseAssembler, dataCollector=iks.io.dataCollectors[1])
    discontinousVtkWriter.addAllResultsAsPointData()

    discontinousVtkWriter.addInterpolation(x, flatBasis, "displacements", 2)
    discontinousVtkWriter.write("file4")



if __name__ == "__main__":
    runTest()
