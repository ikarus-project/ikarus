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


import unittest


class TestVtkWriter(unittest.TestCase):
    def setUp(self):
        lowerLeft = []
        upperRight = []
        elements = []
        for i in range(2):
            lowerLeft.append(-1)
            upperRight.append(1)
            elements.append(3)

        self.grid = dune.grid.structuredGrid(lowerLeft, upperRight, elements)

        basisLagrange1 = iks.basis(
            self.grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
        )
        self.flatBasis = basisLagrange1.flat()
        d = np.zeros(len(self.flatBasis))
        d[0] = 0.0

        lambdaLoad = iks.ValueWrapper(3.0)

        def vL(x, lambdaVal):
            return np.array([lambdaVal * x[0] * 2, 2 * lambdaVal * x[1] * 0])

        vLoad = iks.finite_elements.volumeLoad2D(vL)

        fes = []
        linElastic = iks.finite_elements.linearElastic(youngs_modulus=1000, nu=0.2)
        for e in self.grid.elements:
            fes.append(iks.finite_elements.makeFE(basisLagrange1, linElastic, vLoad))
            fes[-1].bind(e)

        req = fes[0].createRequirement()
        req.insertParameter(lambdaLoad)
        req.insertGlobalSolution(d)

        self.dirichletValues = iks.dirichletValues(self.flatBasis)

        def fixLeftHandEdge(vec, localIndex, localView, intersection):
            if intersection.geometry.center[1] < -0.99999:
                vec[localView.index(localIndex)] = True

        self.dirichletValues.fixBoundaryDOFsUsingLocalViewAndIntersection(
            fixLeftHandEdge
        )

        self.sparseAssembler = iks.assembler.sparseFlatAssembler(
            fes, self.dirichletValues
        )
        self.sparseAssembler.bind(
            req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full
        )

        Msparse = self.sparseAssembler.matrix()
        forces = self.sparseAssembler.vector()

        self.x = sp.sparse.linalg.spsolve(Msparse, -forces)
        req.insertGlobalSolution(self.x)

    def test(self):
        _ = iks.io.vtkWriter(
            self.sparseAssembler,
            datatype=DataTypes.Float64,
            headertype=DataTypes.UInt32,
        )

        writer = iks.io.vtkWriter(self.sparseAssembler, format=FormatTypes.ascii)
        writer.addInterpolationAsPointData(
            self.x, self.flatBasis, "displacements", size=2
        )
        writer.addInterpolationAsCellData(
            self.x, self.flatBasis, "displacements", size=2
        )

        xDisplacementBasis = dune.functions.subspaceBasis(self.flatBasis, 0)
        writer.addInterpolationAsPointData(self.x, xDisplacementBasis, "u", size=1)

        writer.addAllResultsAsCellData()
        writer.addAllResultsAsPointData()

        fileName = writer.write("file")
        self.assertEqual(fileName[-3:], "vtu")

        self.assertEqual(iks.io.dataCollectors[0], "lagrange")
        self.assertEqual(iks.io.dataCollectors[1], "discontinuous")

        writer2 = iks.io.vtkWriter(
            self.sparseAssembler, dataCollector="lagrange", order=2
        )

        writer2.addResultAsCellData("linearStress")
        writer2.addResultAsPointData("linearStress")

        writer2.setFormat(FormatTypes.ascii)
        writer2.setDatatype(DataTypes.Float64)
        writer2.setHeadertype(DataTypes.UInt16)

        @gridFunction(self.grid)
        def g(x):
            return [math.sin(2 * math.pi * x[0] * x[1]), x[0] * x[1]] * 5

        writer2.addPointData(g, name="g", components=(0, 1))
        writer2.addPointData(g, name="g2", components=[0, 1, 2])

        writer2.write("file2")

        # Structured writer
        writer3 = iks.io.vtkWriter(
            self.sparseAssembler, structured=True, format=FormatTypes.ascii
        )
        writer3.addAllResultsAsCellData()
        fileName = writer3.write("file3")
        assert fileName[-3:] == "vtr"

        # These two constructors should result in a warning
        with self.assertWarns(Warning):
            _ = iks.io.vtkWriter(
                self.sparseAssembler, dataCollector="lagrange", structured=True
            )
        with self.assertWarns(Warning):    
            _ = iks.io.vtkWriter(self.sparseAssembler, dataCollector="unknown")

        discontinuousVtkWriter = iks.io.vtkWriter(
            self.sparseAssembler, dataCollector=iks.io.dataCollectors[1]
        )
        discontinuousVtkWriter.addAllResultsAsPointData()

        discontinuousVtkWriter.addInterpolationAsPointData(
            self.x, self.flatBasis, "displacements", size=2
        )
        discontinuousVtkWriter.write("file4")


if __name__ == "__main__":
    unittest.main()
