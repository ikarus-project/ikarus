# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import setpath

setpath.set_path()
import ikarus as iks
import ikarus.finite_elements
import ikarus.utils
import ikarus.assembler
import ikarus.dirichlet_values
import numpy as np
import scipy as sp
from scipy.optimize import minimize

import dune.grid
import dune.functions

if __name__ == "__main__":
    lowerLeft = []
    upperRight = []
    elements = []
    for i in range(2):
        lowerLeft.append(-1)
        upperRight.append(1)
        elements.append(3)

    req = ikarus.FERequirements()
    req.addAffordance(iks.ScalarAffordances.mechanicalPotentialEnergy)

    grid = dune.grid.structuredGrid(lowerLeft, upperRight, elements)
    grid.hierarchicalGrid.globalRefine(0)
    basisLagrange1 = ikarus.basis(
        grid, dune.functions.Power(dune.functions.Lagrange(order=1), 2)
    )
    flatBasis = basisLagrange1.flat()

    forces = np.zeros(8)
    stiffness = np.zeros((8, 8))

    def volumeLoad(x, lambdaVal):
        return np.array([lambdaVal * x[0] * 2, 2 * lambdaVal * x[1] * 0])

    def neumannLoad(x, lambdaVal):
        return np.array([lambdaVal * 0, lambdaVal])

    neumannVertices = np.zeros(grid.size(2) * 2, dtype=bool)
    assert len(neumannVertices) == len(flatBasis)

    flatBasis.interpolate(neumannVertices, lambda x: True if x[1] > 0.9 else False)

    boundaryPatch = iks.utils.boundaryPatch(grid, neumannVertices)

    svk = iks.materials.StVenantKirchhoff(E=1000, nu=0.3)

    psNH = svk.asPlaneStress()
    fes = []
    for e in grid.elements:
        fes.append(
            iks.finite_elements.NonLinearElastic(
                basisLagrange1, e, psNH, volumeLoad, boundaryPatch, neumannLoad
            )
        )

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

    dRed = np.zeros(assembler.reducedSize())

    lambdaLoad = iks.ValueWrapper(3.0)

    def energy(dRedInput):
        reqL = ikarus.FERequirements()
        reqL.addAffordance(iks.ScalarAffordances.mechanicalPotentialEnergy)
        reqL.insertParameter(iks.FEParameter.loadfactor, lambdaLoad)

        dBig = assembler.createFullVector(dRedInput)
        reqL.insertGlobalSolution(iks.FESolutions.displacement, dBig)
        return assembler.getScalar(reqL)

    def gradient(dRedInput):
        reqL = ikarus.FERequirements()
        reqL.addAffordance(iks.VectorAffordances.forces)
        reqL.insertParameter(iks.FEParameter.loadfactor, lambdaLoad)

        dBig = assembler.createFullVector(dRedInput)
        reqL.insertGlobalSolution(iks.FESolutions.displacement, dBig)
        return assembler.getReducedVector(reqL)

    def hess(dRedInput):
        reqL = ikarus.FERequirements()
        reqL.addAffordance(iks.MatrixAffordances.stiffness)
        reqL.insertParameter(iks.FEParameter.loadfactor, lambdaLoad)

        dBig = assembler.createFullVector(dRedInput)
        reqL.insertGlobalSolution(iks.FESolutions.displacement, dBig)
        return assembler.getReducedMatrix(reqL).todense()

    resultd = minimize(energy, x0=dRed, options={"disp": True}, tol=1e-14)
    resultd2 = minimize(
        energy, x0=dRed, jac=gradient, options={"disp": True}, tol=1e-14
    )
    resultd3 = minimize(
        energy,
        method="trust-constr",
        x0=dRed,
        jac=gradient,
        hess=hess,
        options={"disp": True},
    )
    resultd4 = sp.optimize.root(gradient, jac=hess, x0=dRed, tol=1e-10)

    assert np.allclose(resultd.x, resultd2.x, atol=1e-6)
    assert np.allclose(resultd3.x, resultd4.x)
    assert np.all(abs(resultd3.grad) < 1e-8)
    assert np.all(abs(resultd4.fun) < 1e-8)

    req = ikarus.FERequirements()
    fullD = assembler.createFullVector(resultd2.x)
    req.insertGlobalSolution(iks.FESolutions.displacement, fullD)

    res1 = fes[0].calculateAt(req, np.array([0.5, 0.5]), iks.ResultType.PK2Stress)
