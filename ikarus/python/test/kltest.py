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

from dune.iga import ControlPoint, ControlPointNet, NurbsPatchData, IGAGrid
from dune.iga.basis import defaultGlobalBasis, Power, Nurbs

from dune.common.hashit import hashIt
from dune.iga.basis import preBasisTypeName
from dune.generator.generator import SimpleGenerator


def globalBasis(gv, tree):
    generator = SimpleGenerator("Basis", "Ikarus::Python")

    pbfName = preBasisTypeName(tree, gv.cppTypeName)
    element_type = f"Ikarus::Basis<{pbfName}>"
    includes = []
    includes += list(gv.cppIncludes)
    includes += ["dune/iga/nurbsbasis.hh"]
    includes += ["ikarus/python/basis/basis.hh"]

    moduleName = "Basis_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )
    basis = defaultGlobalBasis(gv, tree)
    return module.Basis(basis)


if __name__ == "__main__":
    cp = ControlPoint((0, 0, 0), 1)
    cp2 = ControlPoint((10, 0, 0), 1)
    cp3 = ControlPoint((0, 2, 0), 1)
    cp4 = ControlPoint((10, 2, 0), 1)

    netC = ((cp, cp2), (cp3, cp4))
    net = ControlPointNet(netC)

    nurbsPatchData1 = NurbsPatchData(((0, 0, 1, 1), (0, 0, 1, 1)), net, (1, 1))
    nurbsPatchData2 = nurbsPatchData1.degreeElevate(0, 1)
    nurbsPatchDataFinal = nurbsPatchData2.degreeElevate(1, 1)
    # nurbsPatchDataFinal = nurbsPatchData2.knotRefinement((0.1,0.2,0.3,0.4,0.5), 0)
    # # nurbsPatchData = nurbsPatchData.knotRefinement((0.25,0.5,0.75), 1)

    gridView = IGAGrid(nurbsPatchDataFinal)

    gridView.hierarchicalGrid.globalRefine(2)
    basis = globalBasis(gridView, Power(Nurbs(), 3))

    flatBasis = basis.flat()
    print("Number of dofs: ", len(flatBasis))
    d = np.zeros(len(flatBasis))

    lambdaLoad = iks.ValueWrapper(1.0)
    thickness = 0.1

    def volumeLoad(x, lambdaVal):
        return np.array([0, 0, 2 * thickness**3 * lambdaVal])

    fes = []
    for e in gridView.elements:
        fes.append(
            iks.finite_elements.KirchhoffLoveShell(
                basis, e, 1000, 0.0, thickness, volumeLoad
            )
        )

    dirichletValues = iks.dirichletValues(flatBasis)

    def fixLeftAndRightEdge(vec, localIndex, localView, intersection):
        if (
            intersection.geometry.center[0] < 1e-8
            or intersection.geometry.center[0] > 10 - 1e-8
        ):
            vec[localView.index(localIndex)] = True

    dirichletValues.fixBoundaryDOFsUsingLocalViewAndIntersection(fixLeftAndRightEdge)

    assembler = iks.assembler.sparseFlatAssembler(fes, dirichletValues)

    def gradAndhess(dRedInput):
        req = ikarus.FERequirements()
        req.addAffordance(iks.ScalarAffordances.mechanicalPotentialEnergy)
        req.insertParameter(iks.FEParameter.loadfactor, lambdaLoad)
        dBig = assembler.createFullVector(dRedInput)
        req.insertGlobalSolution(iks.FESolutions.displacement, dBig)
        g = assembler.getReducedVector(req)
        h = assembler.getReducedMatrix(req)
        return [g, h]

    from numpy.linalg import norm

    maxiter = 100
    abs_tolerance = 1e-8
    d = np.zeros(assembler.reducedSize())
    for k in range(maxiter):
        R, K = gradAndhess(d)
        r_norm = norm(R)

        deltad = sp.sparse.linalg.spsolve(K, R)
        d -= deltad
        print(k, r_norm, norm(deltad))
        if r_norm < abs_tolerance:
            break
    print(d)
    dBig = assembler.createFullVector(d)
    displacementFunc = flatBasis.asFunction(dBig)
    vtkWriter = gridView.trimmedVtkWriter()
    vtkWriter.addPointData(displacementFunc, name="displacement")
    vtkWriter.write("KLshell")
    assert (0.2087577577980777 - max(d)) < 1e-6, f"The maximum displacement should be 0.2087577577980777 but is {max(d)}"
