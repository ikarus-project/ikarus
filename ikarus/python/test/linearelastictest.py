# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import debug_info

debug_info.setDebugFlags()

import ikarus as iks
from ikarus import finite_elements, utils, assembler
import numpy as np
import scipy as sp

import dune.grid
import dune.functions
from dune.vtk import vtkWriter, vtkUnstructuredGridWriter


def linElasticTest(fe_type="standard"):
    """
    fe_type: "standard" (default), "eas" or "assumed"
        When fe_type is "eas" or "assumed", the test is in mixed mode (isMixed=True)
        and the corresponding finite element (eas(4) or assumedStress(5)) is used.
    """
    # Determine mixed formulation option
    isMixed = fe_type in ("eas", "assumed")
    feOption = None
    if fe_type == "eas":
        feOption = iks.finite_elements.eas(4)
    elif fe_type == "assumed":
        feOption = iks.finite_elements.assumedStress(5)

    # Create grid and bases
    lowerLeft = [-1, -1]
    upperRight = [1, 1]
    elements = [3, 3]
    grid = dune.grid.structuredGrid(lowerLeft, upperRight, elements)

    basisLagrange1 = iks.basis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    )
    flatBasis = basisLagrange1.flat()
    d = np.zeros(len(flatBasis))
    d[0] = 0.0

    # Define loads and load functions
    lambdaLoad = iks.Scalar(3.0)

    def vL(x, lambdaVal):
        return np.array([lambdaVal * x[0] * 2, 0])

    vLoad = iks.finite_elements.volumeLoad2D(vL)

    def neumannLoad(x, lambdaVal):
        return np.array([0, lambdaVal])

    neumannVertices = np.zeros(grid.size(2), dtype=bool)

    def loadTopEdgePredicate(x):
        return x[1] > 0.999

    indexSet = grid.indexSet
    for v in grid.vertices:
        neumannVertices[indexSet.index(v)] = loadTopEdgePredicate(v.geometry.center)

    boundaryPatch = iks.utils.boundaryPatch(grid, neumannVertices)
    nBLoad = iks.finite_elements.neumannBoundaryLoad(boundaryPatch, neumannLoad)

    linMat = iks.materials.LinearElasticity(E=1000, nu=0.2).asPlaneStress()
    linElastic = iks.finite_elements.linearElastic(linMat)

    # Assemble the finite elements
    fes = []
    for e in grid.elements:
        if isMixed:
            fe = iks.finite_elements.makeFE(
                basisLagrange1, linElastic, feOption, vLoad, nBLoad
            )
        else:
            fe = iks.finite_elements.makeFE(basisLagrange1, linElastic, vLoad, nBLoad)
        fe.bind(e)
        fes.append(fe)

    # Create requirement object, insert parameters and global solution vectors
    req = fes[0].createRequirement()
    req.insertParameter(lambdaLoad)
    assert req.parameter() == lambdaLoad
    req.insertGlobalSolution(d)

    # Calculate forces and stiffness for first FE and get its local view
    forces = np.zeros(8)
    stiffness = np.zeros((8, 8))
    fes[0].calculateVector(req, iks.VectorAffordance.forces, forces)
    fes[0].calculateMatrix(req, iks.MatrixAffordance.stiffness, stiffness)
    fes[0].localView()

    # Create additional Dirichlet values and fix boundary DOFs using helper callbacks
    dirichletValues = iks.dirichletValues(flatBasis)

    def fixFirstIndex(vec, globalIndex):
        vec[0] = True

    def fixAnotherVertex(vec, localIndex, localView):
        localView.index(localIndex)
        vec[1] = True

    def fixLeftHandEdge(vec, localIndex, localView, intersection):
        if intersection.geometry.center[1] < -0.99999:
            vec[localView.index(localIndex)] = True

    dirichletValues.fixBoundaryDOFs(fixFirstIndex)
    dirichletValues.fixBoundaryDOFs(fixAnotherVertex)
    dirichletValues.fixBoundaryDOFs(fixLeftHandEdge)

    # Assemble matrices and vectors
    assembler = iks.assembler.sparseFlatAssembler(fes, dirichletValues)
    assembler.bind(req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full)

    Msparse = assembler.matrix()
    forces = assembler.vector()

    x = sp.sparse.linalg.spsolve(Msparse, -forces)
    fx = flatBasis.asFunction(x)
    req.globalSolution()

    # Setup stress functions for VTK writing
    stressFuncScalar = grid.function(
        lambda e, x: fes[indexSet.index(e)].calculateAt(req, x, "linearStress")[0]
    )
    stressFuncVec = grid.function(
        lambda e, x: fes[indexSet.index(e)].calculateAt(req, x, "linearStress")[:]
    )

    fileName = "resultdisplacement" + (
        "EAS" if fe_type == "eas" else ("AssumedStress" if fe_type == "assumed" else "")
    )
    vtkWriter(grid, fileName, pointData={("displacement", (0, 1)): fx})

    writer2 = vtkUnstructuredGridWriter(grid)
    writer2.addCellData(stressFuncScalar, name="stress")
    writer2.addCellData(stressFuncVec, name="stress2")
    writer2.write(
        name="result"
        + (
            "EAS"
            if fe_type == "eas"
            else ("AssumedStress" if fe_type == "assumed" else "")
        )
    )

    # Querying for an unsupported ResultType should trigger a runtime error.
    try:
        fes[0].calculateAt(req, np.array([0.5, 0.5]), "PK2Stress")
    except RuntimeError:
        pass
    else:
        raise AssertionError("Expected RuntimeError for PK2Stress query")

    # We only test assembler manipulators here for standard element, as there is nothing specific about the combination of
    # element technology and manipulators and it only adds compile time
    if not isMixed:
        assemblerDense = iks.assembler.denseFlatAssembler(fes, dirichletValues)
        assemblerDense.bind(
            req, iks.AffordanceCollection.elastoStatics, iks.DBCOption.Full
        )

        # Create assembler manipulators and add call-backs.
        assemblerManipulator = iks.assembler.assemblerManipulator(assembler)
        assemblerManipulatorDense = iks.assembler.assemblerManipulator(assemblerDense)

        def scalarf(assembler, req, affordance, scalar):
            scalar *= 2

        assemblerManipulator.addScalarCallBack(scalarf)
        assemblerManipulatorDense.addScalarCallBack(scalarf)
        req3 = assemblerManipulator.requirement()
        d3 = req3.globalSolution()
        d3 += 1
        assert (
            abs(2 * assembler.scalar() - assemblerManipulator.scalar()) < 1e-6
        ), f"assembler.scalar(): {assembler.scalar()}, assemblerManipulator.scalar(): {assemblerManipulator.scalar()}"
        assert (
            abs(2 * assembler.scalar() - assemblerManipulatorDense.scalar()) < 1e-6
        ), f"assembler.scalar(): {assembler.scalar()}, assemblerManipulator.scalar(): {assemblerManipulatorDense.scalar()}"

        def vectorf(assembler, req, affordance, dbcOption, vector):
            vector[5] += 2

        assemblerManipulator.addVectorCallBack(vectorf)
        assemblerManipulatorDense.addVectorCallBack(vectorf)
        assert (
            abs(assembler.vector()[5] - (assemblerManipulator.vector()[5] - 2)) < 1e-6
        ), f"assembler.vector()[5]: {assembler.vector()[5]}, assemblerManipulator.vector()[5]: {assemblerManipulator.vector()[5]}"
        assert (
            abs(assembler.vector()[5] - (assemblerManipulatorDense.vector()[5] - 2))
            < 1e-6
        ), f"assembler.vector()[5]: {assembler.vector()[5]}, assemblerManipulator.vector()[5]: {assemblerManipulatorDense.vector()[5]}"

        def matrixf(assembler, req, affordance, dbcOption, matrix):
            matrix[5, 6] += 2.0

        assemblerManipulator.addMatrixCallBack(matrixf)
        assemblerManipulatorDense.addMatrixCallBack(matrixf)
        assert (
            abs(assembler.matrix()[5, 6] - (assemblerManipulator.matrix()[5, 6] - 2))
            < 1e-6
        ), f"assembler.matrix()[5,6]: {assembler.matrix()[5,6]}, assemblerManipulator.matrix()[5,6]: {assemblerManipulator.matrix()[5,6]}"
        assert (
            abs(
                assembler.matrix()[5, 6]
                - (assemblerManipulatorDense.matrix()[5, 6] - 2)
            )
            < 1e-6
        ), f"assembler.matrix()[5,6]: {assembler.matrix()[5,6]}, assemblerManipulator.matrix()[5,6]: {assemblerManipulatorDense.matrix()[5,6]}"

        # Test global index mapping on the structured grid
        startX = -1
        incX = 2 / 3
        startY = -1
        incY = 2 / 3
        for j in range(4):
            for i in range(4):
                indices = iks.utils.globalIndexFromGlobalPosition(
                    flatBasis, (startX + i * incX, startY + j * incY)
                )
                assert (
                    indices[0] == i + 4 * j
                ), f"Expected {i + 4 * j} from {i},{j}, got {indices[0]} at pos {(startX+i*incX, startY+j*incY)}"
                assert (
                    indices[1] == i + 4 * j + 16
                ), f"Expected {i + 4 * j + 16} from {i},{j}, got {indices[1]} at pos {(startX+i*incX, startY+j*incY)}"

        # Final assertions on assembler sizes and constraints
        assert assembler.size == assemblerManipulator.size == flatBasis.size()
        assert len(assembler) == len(assemblerManipulator) == len(flatBasis)
        assert (
            assembler.isConstrained(1)
            == assemblerManipulator.isConstrained(1)
            == dirichletValues.isConstrained(1)
        )
        # Ensure constraintsBelow returns matching values
        for i in range(flatBasis.size()):
            assert assembler.constraintsBelow(
                i
            ) == assemblerManipulator.constraintsBelow(
                i
            ), f"Mismatch in constraintsBelow for index {i}"


if __name__ == "__main__":
    # # Run for standard formulation.
    linElasticTest(fe_type="standard")

    # Run for eas formulation.
    linElasticTest(fe_type="eas")

    # Run for assumed stress formulation.
    linElasticTest(fe_type="assumed")
