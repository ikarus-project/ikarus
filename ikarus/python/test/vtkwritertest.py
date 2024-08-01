# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# import debug_info
# debug_info.setDebugFlags()

import math

import ikarus as iks
from ikarus import finite_elements, utils, assembler, io

import numpy as np
import scipy as sp

import dune.grid
import dune.functions

from dune.grid import gridFunction
from dune.vtk import FormatTypes, DataTypes

import os
import unittest

from dirichletvaluetest import makeGrid


class TestVtkWriter(unittest.TestCase):
    def setUp(self):
        reader = (
            dune.grid.reader.gmsh,
            os.path.join(os.path.dirname(__file__), "auxiliaryfiles/quad2d.msh"),
        )

        self.grid = dune.grid.ugGrid(reader, dimgrid=2)

        basis = iks.basis(
            self.grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
        )
        self.flatBasis = basis.flat()
        d = np.zeros(len(self.flatBasis))

        lambdaLoad = iks.Scalar(3.0)

        def vL(x, lambdaVal):
            return np.array([lambdaVal * x[0] * 2, 2 * lambdaVal * x[1] * 0])

        vLoad = iks.finite_elements.volumeLoad2D(vL)

        fes = []
        linElastic = iks.finite_elements.linearElastic(youngs_modulus=1000, nu=0.2)
        for e in self.grid.elements:
            fes.append(iks.finite_elements.makeFE(basis, linElastic, vLoad))
            fes[-1].bind(e)

        req = fes[0].createRequirement()
        req.insertParameter(lambdaLoad)
        req.insertGlobalSolution(d)

        self.dirichletValues = iks.dirichletValues(self.flatBasis)

        def fixLeftHandEdge(vec, localIndex, localView, intersection):
            if intersection.geometry.center[1] < -0.99999:
                vec[localView.index(localIndex)] = True

        self.dirichletValues.fixBoundaryDOFs(fixLeftHandEdge)

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
        writer.addInterpolation(
            self.x, self.flatBasis, "displacements", io.DataTag.asPointData
        )
        writer.addInterpolation(
            self.x, self.flatBasis, "displacements", io.DataTag.asCellData
        )

        xDisplacementBasis = dune.functions.subspaceBasis(self.flatBasis, 0)
        writer.addInterpolation(self.x, xDisplacementBasis, "u", io.DataTag.asPointData)

        writer.addAllResults(io.DataTag.asCellData)
        writer.addAllResults(io.DataTag.asPointData)

        fileName = writer.write("file")
        self.assertEqual(fileName[-3:], "vtu")

    def test_lagrange(self):
        writer2 = iks.io.vtkWriter(
            self.sparseAssembler, dataCollector=iks.io.DataCollector.lagrange, order=2
        )

        writer2.addResult("linearStress", io.DataTag.asCellData)
        writer2.addResult("linearStress", io.DataTag.asPointData)

        writer2.setFormat(FormatTypes.ascii)
        writer2.setDatatype(DataTypes.Float64)
        writer2.setHeadertype(DataTypes.UInt16)

        @gridFunction(self.grid)
        def g(x):
            return [math.sin(2 * math.pi * x[0] * x[1]), x[0] * x[1]] * 5

        writer2.addPointData(g, name="g", components=(0, 1))
        writer2.addPointData(g, name="g2", components=[0, 1, 2])

        writer2.write("file2")

    def test_structuredWriter(self):
        gridUG = makeGrid()

        basis = iks.basis(
            gridUG, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
        )
        flatBasisYASP = basis.flat()
        d = np.zeros(len(flatBasisYASP))

        lambdaLoad = iks.Scalar(3.0)

        fes = []
        linElastic = iks.finite_elements.linearElastic(youngs_modulus=1000, nu=0.2)
        for e in gridUG.elements:
            fes.append(iks.finite_elements.makeFE(basis, linElastic))
            fes[-1].bind(e)

        req = fes[0].createRequirement()
        lambdaLoad = iks.Scalar(3.0)
        req.insertParameter(lambdaLoad)
        req.insertGlobalSolution(d)

        dirichletValuesYASP = iks.dirichletValues(flatBasisYASP)
        sparseAssemblerYASP = iks.assembler.sparseFlatAssembler(
            fes, dirichletValuesYASP
        )

        sparseAssemblerYASP.bind(
            req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full
        )

        writer3 = iks.io.vtkWriter(sparseAssemblerYASP, format=FormatTypes.ascii)
        writer3.addAllResults(io.DataTag.asCellData)

        print(writer3.cppTypeName)

        fileName = writer3.write("file3")
        self.assertEqual(fileName[-3:], "vtr")

        # with self.assertWarns(Warning):
        #     _ = iks.io.vtkWriter(self.sparseAssembler, dataCollector="unknown")

    def test_discontinuousWriter(self):
        discontinuousVtkWriter = iks.io.vtkWriter(
            self.sparseAssembler, dataCollector=iks.io.DataCollector.discontinuous
        )
        discontinuousVtkWriter.addAllResults(io.DataTag.asPointData)

        discontinuousVtkWriter.addInterpolation(
            self.x, self.flatBasis, "displacements", io.DataTag.asPointData
        )
        discontinuousVtkWriter.write("file4")


if __name__ == "__main__":
    unittest.main()
