# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import setpath

setpath.set_path()
import ikarus as iks
import numpy as np
import math


def checkSizes(A, expected_rows, expected_cols):
    shape = A.shape
    if len(shape) == 2:
        num_rows, num_cols = shape
    else:
        num_rows = shape[0]
        num_cols = 1

    assert expected_rows == num_rows
    assert expected_cols == num_cols


def checkMaterial(mat, strain, nonlinear=True, size=6):
    if nonlinear:
        mat.storedEnergy("greenLagrangian", strain)
        mat.storedEnergy("rightCauchyGreenTensor", strain)
        S = mat.stresses("greenLagrangian", strain)
        checkSizes(S, size, 1)
        S = mat.stresses("rightCauchyGreenTensor", strain)
        checkSizes(S, size, 1)
        C = mat.tangentModuli("greenLagrangian", strain)
        checkSizes(C, size, size)
        C = mat.tangentModuli("rightCauchyGreenTensor", strain)
        checkSizes(C, size, size)
        C = mat.tangentModuli("greenLagrangian", strain)
        checkSizes(C, size, size)
    else:
        mat.storedEnergy("linear", strain)
        S = mat.stresses("linear", strain)
        checkSizes(S, size, 1)
        C = mat.tangentModuli("linear", strain)
        checkSizes(C, size, size)

    try:
        mat.storedEnergy(strain)
        assert False
    except TypeError:
        pass

    try:
        mat.stresses("displacementGradient", strain)
        assert False
    except RuntimeError:
        pass

    try:
        mat.stresses("deformationGradient", strain)
        assert False
    except RuntimeError:
        pass


def checkPlaneStressReducedFullEquality(material, strainName, strain):
    strainFull = np.array(
        [
            strain[0],
            strain[1],
            0 if strainName == "linear" or strainName == "greenLagrangian" else 1,
            0,
            0,
            strain[2],
        ]
    )

    e1 = material.storedEnergy(strainName, strainFull)
    e2 = material.storedEnergy(strainName, strain)

    assert math.isclose(e1, e2)

    s1 = material.stresses(strainName, strainFull)
    s2 = material.stresses(strainName, strain)
    assert np.allclose(s1, s2)

    C1 = material.tangentModuli(strainName, strainFull)
    C2 = material.tangentModuli(strainName, strain)
    assert np.allclose(C1, C2)


def checkShellStressReducedFullEquality(material, strainName, strain):
    strainFull = np.array(
        [
            strain[0],
            strain[1],
            0 if strainName == "linear" or strainName == "greenLagrangian" else 1,
            strain[2],
            strain[3],
            strain[4],
        ]
    )

    e1 = material.storedEnergy(strainName, strainFull)
    e2 = material.storedEnergy(strainName, strain)

    assert math.isclose(e1, e2)

    s1 = material.stresses(strainName, strainFull)
    s2 = material.stresses(strainName, strain)
    assert np.allclose(s1, s2)

    C1 = material.tangentModuli(strainName, strainFull)
    C2 = material.tangentModuli(strainName, strain)
    assert np.allclose(C1, C2)


def checkBeamStressReducedFullEquality(material, strainName, strain):
    strainFull = np.array(
        [
            strain[0],
            0 if strainName == "linear" or strainName == "greenLagrangian" else 1,
            0 if strainName == "linear" or strainName == "greenLagrangian" else 1,
            strain[1],
            strain[2],
            strain[3],
        ]
    )

    e1 = material.storedEnergy(strainName, strainFull)
    e2 = material.storedEnergy(strainName, strain)

    assert math.isclose(e1, e2)

    s1 = material.stresses(strainName, strainFull)
    s2 = material.stresses(strainName, strain)
    assert np.allclose(s1, s2)

    C1 = material.tangentModuli(strainName, strainFull)
    C2 = material.tangentModuli(strainName, strain)
    assert np.allclose(C1, C2)


def checkWithStrain(strain):
    svk = iks.materials.StVenantKirchhoff(E=1000, nu=0.3)
    nh = iks.materials.NeoHooke(E=1000, nu=0.3)
    lin = iks.materials.LinearElasticity(E=1000, nu=0.3)

    if len(strain) == 6:
        checkMaterial(svk, strain)
        checkMaterial(nh, strain)
        checkMaterial(lin, strain, False)
        checkMaterial(lin.asPlaneStress(), strain, False, 3)
        checkMaterial(nh.asPlaneStress(), strain, True, 3)
        checkMaterial(svk.asPlaneStress(), strain, True, 3)

        checkMaterial(lin.asShellMaterial(), strain, False, 5)
        checkMaterial(nh.asShellMaterial(), strain, True, 5)
        checkMaterial(svk.asShellMaterial(), strain, True, 5)

        checkMaterial(lin.asBeamMaterial(), strain, False, 4)
        checkMaterial(nh.asBeamMaterial(), strain, True, 4)
        checkMaterial(svk.asBeamMaterial(), strain, True, 4)
    elif len(strain) == 3:
        checkMaterial(lin.asPlaneStress(), strain, False, 3)
        checkMaterial(nh.asPlaneStress(), strain, True, 3)
        checkMaterial(svk.asPlaneStress(), strain, True, 3)

        checkPlaneStressReducedFullEquality(
            nh.asPlaneStress(), "rightCauchyGreenTensor", strain
        )
        checkPlaneStressReducedFullEquality(
            svk.asPlaneStress(), "rightCauchyGreenTensor", strain
        )

        strain[0] = strain[0] - 1
        strain[1] = strain[1] - 1
        checkPlaneStressReducedFullEquality(lin.asPlaneStress(), "linear", strain)
        checkPlaneStressReducedFullEquality(
            nh.asPlaneStress(), "greenLagrangian", strain
        )
        checkPlaneStressReducedFullEquality(
            svk.asPlaneStress(), "greenLagrangian", strain
        )

    elif len(strain) == 5:
        checkMaterial(lin.asShellMaterial(), strain, False, 5)
        checkMaterial(nh.asShellMaterial(), strain, True, 5)
        checkMaterial(svk.asShellMaterial(), strain, True, 5)

        checkShellStressReducedFullEquality(
            nh.asShellMaterial(), "rightCauchyGreenTensor", strain
        )
        checkShellStressReducedFullEquality(
            svk.asShellMaterial(), "rightCauchyGreenTensor", strain
        )

        strain[0] = strain[0] - 1
        strain[1] = strain[1] - 1
        checkShellStressReducedFullEquality(lin.asShellMaterial(), "linear", strain)
        checkShellStressReducedFullEquality(
            nh.asShellMaterial(), "greenLagrangian", strain
        )
        checkShellStressReducedFullEquality(
            svk.asShellMaterial(), "greenLagrangian", strain
        )

    elif len(strain) == 4:
        checkMaterial(lin.asBeamMaterial(), strain, False, 4)
        checkMaterial(nh.asBeamMaterial(), strain, True, 4)
        checkMaterial(svk.asBeamMaterial(), strain, True, 4)

        checkBeamStressReducedFullEquality(
            nh.asBeamMaterial(), "rightCauchyGreenTensor", strain
        )
        checkBeamStressReducedFullEquality(
            svk.asBeamMaterial(), "rightCauchyGreenTensor", strain
        )

        strain[0] = strain[0] - 1
        checkBeamStressReducedFullEquality(lin.asBeamMaterial(), "linear", strain)
        checkBeamStressReducedFullEquality(
            nh.asBeamMaterial(), "greenLagrangian", strain
        )
        checkBeamStressReducedFullEquality(
            svk.asBeamMaterial(), "greenLagrangian", strain
        )


if __name__ == "__main__":
    help(iks)
    help(iks.materials)

    strain = np.array([1.2, 1.1, 0.9, 0.1, 0.2, 0.2])
    checkWithStrain(strain)

    # check if passing plane strains is enough
    strain = np.array([1.2, 1.1, 0.1])
    checkWithStrain(strain)

    # check if passing shell strains is enough
    strain = np.array([1.2, 1.1, 0.1, 0.05, 0.2])
    checkWithStrain(strain)

    # check if passing beam strains is enough
    strain = np.array([1.2, 0.1, 0.05, 0.2])
    checkWithStrain(strain)
