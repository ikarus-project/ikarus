# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import debug_info

debug_info.setDebugFlags()

import ikarus as iks
from ikarus import materials, utils

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


def checkMaterial(mat, strain, strainTag, size=6):
    mat.storedEnergy(strainTag, strain)
    S = mat.stresses(strainTag, strain)
    checkSizes(S, size, 1)
    C = mat.tangentModuli(strainTag, strain)
    checkSizes(C, size, size)

    try:
        mat.storedEnergy(strain)
        assert False
    except TypeError:
        pass

    try:
        mat.stresses(materials.StrainTags.displacementGradient, strain)
        assert False
    except RuntimeError:
        pass

    try:
        mat.stresses(materials.StrainTags.deformationGradient, strain)
        assert False
    except RuntimeError:
        pass


# This checks plane stress and plane strain
def check2DReducedFullEquality(material, strainName, strain):
    strainFull = np.array(
        [
            strain[0],
            strain[1],
            (
                0
                if strainName == materials.StrainTags.linear
                or strainName == materials.StrainTags.greenLagrangian
                else 1
            ),
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
            (
                0
                if strainName == materials.StrainTags.linear
                or strainName == materials.StrainTags.greenLagrangian
                else 1
            ),
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
            (
                0
                if strainName == materials.StrainTags.linear
                or strainName == materials.StrainTags.greenLagrangian
                else 1
            ),
            (
                0
                if strainName == materials.StrainTags.linear
                or strainName == materials.StrainTags.greenLagrangian
                else 1
            ),
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
    mlin = materials.muesliMaterial(
        materials.MuesliSmallStrain.elasticIsotropicMaterial, E=100, nu=0.3
    )
    mnh = materials.muesliMaterial(
        materials.MuesliFiniteStrain.neohookeanMaterial, E=1000, nu=0.3
    )
    mmoon = materials.muesliMaterial(
        materials.MuesliFiniteStrain.mooneyMaterial, alpha=[500, 200, 200]
    )
    rCGT = materials.StrainTags.rightCauchyGreenTensor
    linT = materials.StrainTags.linear
    gtT = materials.StrainTags.greenLagrangian

    if len(strain) == 6:
        import copy

        strainS = copy.deepcopy(strain)
        strainS[0] = strain[0] - 1
        strainS[1] = strain[1] - 1
        strainS[2] = strain[2] - 1
        for strainTag in [rCGT, gtT]:
            strainToUse = strain if strainTag is rCGT else strainS
            checkMaterial(svk, strainToUse, strainTag)
            checkMaterial(nh, strainToUse, strainTag)
            checkMaterial(mnh, strainToUse, strainTag)
            checkMaterial(mmoon, strainToUse, strainTag)

            checkMaterial(nh.asPlaneStress(), strainToUse, strainTag, 3)
            checkMaterial(svk.asPlaneStress(), strainToUse, strainTag, 3)
            checkMaterial(mnh.asPlaneStress(), strainToUse, strainTag, 3)
            checkMaterial(mmoon.asPlaneStress(), strainToUse, strainTag, 3)

            checkMaterial(nh.asPlaneStrain(), strainToUse, strainTag, 3)
            checkMaterial(svk.asPlaneStrain(), strainToUse, strainTag, 3)
            checkMaterial(mnh.asPlaneStrain(), strainToUse, strainTag, 3)
            checkMaterial(mmoon.asPlaneStrain(), strainToUse, strainTag, 3)

            checkMaterial(nh.asShellMaterial(), strainToUse, strainTag, 5)
            checkMaterial(svk.asShellMaterial(), strainToUse, strainTag, 5)
            checkMaterial(mnh.asShellMaterial(), strainToUse, strainTag, 5)
            checkMaterial(mmoon.asShellMaterial(), strainToUse, strainTag, 5)

            checkMaterial(nh.asBeamMaterial(), strainToUse, strainTag, 4)
            checkMaterial(svk.asBeamMaterial(), strainToUse, strainTag, 4)
            checkMaterial(mnh.asBeamMaterial(), strainToUse, strainTag, 4)
            checkMaterial(mmoon.asBeamMaterial(), strainToUse, strainTag, 4)

        checkMaterial(lin, strain, linT)
        checkMaterial(mlin, strain, linT)

        checkMaterial(lin.asPlaneStress(), strain, linT, 3)
        checkMaterial(mlin.asPlaneStress(), strain, linT, 3)

        checkMaterial(lin.asPlaneStrain(), strain, linT, 3)
        checkMaterial(mlin.asPlaneStrain(), strain, linT, 3)

        checkMaterial(lin.asShellMaterial(), strain, linT, 5)
        checkMaterial(mlin.asShellMaterial(), strain, linT, 5)

        checkMaterial(lin.asBeamMaterial(), strain, linT, 4)
        checkMaterial(mlin.asBeamMaterial(), strain, linT, 4)

    elif len(strain) == 3:
        for strainTag in [rCGT, gtT]:
            checkMaterial(nh.asPlaneStress(), strain, strainTag, 3)
            checkMaterial(svk.asPlaneStress(), strain, strainTag, 3)

            checkMaterial(nh.asPlaneStrain(), strain, strainTag, 3)
            checkMaterial(svk.asPlaneStrain(), strain, strainTag, 3)

            checkMaterial(nh.asPlaneStrain(), strain, strainTag, 3)
            checkMaterial(svk.asPlaneStrain(), strain, strainTag, 3)

        checkMaterial(mnh.asPlaneStress(), strain, rCGT, 3)
        checkMaterial(mmoon.asPlaneStress(), strain, rCGT, 3)

        checkMaterial(lin.asPlaneStress(), strain, linT, 3)
        checkMaterial(mlin.asPlaneStress(), strain, linT, 3)

        checkMaterial(lin.asPlaneStrain(), strain, linT, 3)
        checkMaterial(mlin.asPlaneStrain(), strain, linT, 3)

        check2DReducedFullEquality(
            nh.asPlaneStress(), materials.StrainTags.rightCauchyGreenTensor, strain
        )
        check2DReducedFullEquality(
            nh.asPlaneStrain(), materials.StrainTags.rightCauchyGreenTensor, strain
        )
        check2DReducedFullEquality(
            svk.asPlaneStress(), materials.StrainTags.rightCauchyGreenTensor, strain
        )
        check2DReducedFullEquality(
            svk.asPlaneStrain(), materials.StrainTags.rightCauchyGreenTensor, strain
        )

        strain[0] = strain[0] - 1
        strain[1] = strain[1] - 1
        check2DReducedFullEquality(
            lin.asPlaneStress(), materials.StrainTags.linear, strain
        )
        check2DReducedFullEquality(
            lin.asPlaneStrain(), materials.StrainTags.linear, strain
        )

        check2DReducedFullEquality(
            nh.asPlaneStress(), materials.StrainTags.greenLagrangian, strain
        )
        check2DReducedFullEquality(
            nh.asPlaneStrain(), materials.StrainTags.greenLagrangian, strain
        )
        check2DReducedFullEquality(
            svk.asPlaneStress(), materials.StrainTags.greenLagrangian, strain
        )
        check2DReducedFullEquality(
            svk.asPlaneStrain(), materials.StrainTags.greenLagrangian, strain
        )

    elif len(strain) == 5:
        for strainTag in [rCGT, gtT]:
            checkMaterial(nh.asShellMaterial(), strain, strainTag, 5)
            checkMaterial(svk.asShellMaterial(), strain, strainTag, 5)

        checkMaterial(lin.asShellMaterial(), strain, linT, 5)
        checkMaterial(mlin.asShellMaterial(), strain, linT, 5)

        checkMaterial(mnh.asShellMaterial(), strain, rCGT, 5)
        checkMaterial(mmoon.asShellMaterial(), strain, rCGT, 5)

        checkShellStressReducedFullEquality(
            nh.asShellMaterial(), materials.StrainTags.rightCauchyGreenTensor, strain
        )
        checkShellStressReducedFullEquality(
            svk.asShellMaterial(), materials.StrainTags.rightCauchyGreenTensor, strain
        )

        strain[0] = strain[0] - 1
        strain[1] = strain[1] - 1
        checkShellStressReducedFullEquality(
            lin.asShellMaterial(), materials.StrainTags.linear, strain
        )
        checkShellStressReducedFullEquality(
            nh.asShellMaterial(), materials.StrainTags.greenLagrangian, strain
        )
        checkShellStressReducedFullEquality(
            svk.asShellMaterial(), materials.StrainTags.greenLagrangian, strain
        )

    elif len(strain) == 4:
        for strainTag in [rCGT, gtT]:
            checkMaterial(nh.asBeamMaterial(), strain, strainTag, 4)
            checkMaterial(svk.asBeamMaterial(), strain, strainTag, 4)

        checkMaterial(lin.asBeamMaterial(), strain, linT, 4)
        checkMaterial(mlin.asBeamMaterial(), strain, linT, 4)

        checkMaterial(mnh.asBeamMaterial(), strain, rCGT, 4)
        checkMaterial(mmoon.asBeamMaterial(), strain, rCGT, 4)

        checkBeamStressReducedFullEquality(
            nh.asBeamMaterial(), materials.StrainTags.rightCauchyGreenTensor, strain
        )
        checkBeamStressReducedFullEquality(
            svk.asBeamMaterial(), materials.StrainTags.rightCauchyGreenTensor, strain
        )

        strain[0] = strain[0] - 1
        checkBeamStressReducedFullEquality(
            lin.asBeamMaterial(), materials.StrainTags.linear, strain
        )
        checkBeamStressReducedFullEquality(
            nh.asBeamMaterial(), materials.StrainTags.greenLagrangian, strain
        )
        checkBeamStressReducedFullEquality(
            svk.asBeamMaterial(), materials.StrainTags.greenLagrangian, strain
        )


def checkStrainTransformation():
    glVoigt = [1.2, 1.1, 0.9, 0.1, 0.2, 0.2]
    glTensor = utils.fromVoigt(glVoigt)

    c = materials.transformStrain(
        materials.StrainTags.greenLagrangian,
        materials.StrainTags.rightCauchyGreenTensor,
        glVoigt,
    )
    c2 = materials.transformStrain(
        materials.StrainTags.greenLagrangian,
        materials.StrainTags.rightCauchyGreenTensor,
        glTensor,
    )

    assert np.allclose(c, 2 * glTensor + np.identity(3))
    assert np.allclose(c2, 2 * glTensor + np.identity(3))

    gl = materials.transformStrain(
        materials.StrainTags.rightCauchyGreenTensor,
        materials.StrainTags.greenLagrangian,
        c,
    )

    assert np.allclose(gl, glTensor)
    assert np.allclose(gl, 0.5 * (c - np.identity(3)))

    gl2 = materials.transformStrain(
        materials.StrainTags.greenLagrangian, materials.StrainTags.greenLagrangian, gl
    )

    assert np.allclose(gl, gl2)


def checkVoigtTransformations():
    # Check Voigt transformations for 2D
    glVoigt = [1.2, 1.1, 0.1]
    glTensor = utils.fromVoigt(glVoigt)

    assert np.allclose(utils.toVoigt(glTensor), glVoigt)
    assert np.allclose(utils.fromVoigt(glVoigt), glTensor)

    # Check Voigt transformations for 1D
    glVoigt = [1.2]
    glTensor = utils.fromVoigt(glVoigt)

    assert np.allclose(utils.toVoigt(glTensor), glVoigt)
    assert np.allclose(utils.fromVoigt(glVoigt), glTensor)

    # Check Voigt transformations for 3D
    glVoigt = [1.2, 1.1, 0.9, 0.1, 0.2, 0.2]
    glTensor = utils.fromVoigt(glVoigt)

    assert np.allclose(utils.toVoigt(glTensor), glVoigt)
    assert np.allclose(utils.fromVoigt(glVoigt), glTensor)


def checkMaterialConstructors():
    # Check different constructors (no physical meaning)
    iks.materials.StVenantKirchhoff(E=1000, mu=500)
    iks.materials.StVenantKirchhoff(E=1000, nu=0.3)
    iks.materials.StVenantKirchhoff(E=1000, K=500)
    iks.materials.StVenantKirchhoff(E=1000, Lambda=500)
    iks.materials.StVenantKirchhoff(K=1000, Lambda=500)
    iks.materials.StVenantKirchhoff(Lambda=1000, mu=500)

    svk1 = iks.materials.StVenantKirchhoff(E=1000, nu=0.3)
    svk2 = iks.materials.StVenantKirchhoff(nu=0.3, E=1000)

    strain = np.array([1.2, 1.1, 0.9, 0.1, 0.2, 0.2])
    s1 = svk1.stresses(materials.StrainTags.greenLagrangian, strain)
    s2 = svk2.stresses(materials.StrainTags.greenLagrangian, strain)

    assert np.allclose(s1, s2)


def instantiateMuesliMaterials():

    materials.muesliMaterial(
        materials.MuesliSmallStrain.elasticIsotropicMaterial, E=100, nu=0.3
    )
    materials.muesliMaterial(materials.MuesliFiniteStrain.svkMaterial, E=100, nu=0.3)
    materials.muesliMaterial(
        materials.MuesliFiniteStrain.neohookeanMaterial, E=100, nu=0.3
    )
    materials.muesliMaterial(
        materials.MuesliFiniteStrain.mooneyMaterial, alpha=[500, 200, 200]
    )
    materials.muesliMaterial(
        materials.MuesliFiniteStrain.mooneyMaterial,
        alpha=[0.2, 0.3, 0.4],
        incompressible=True,
    )
    materials.muesliMaterial(
        materials.MuesliFiniteStrain.arrudaboyceMaterial, C1=500, lambda_m=200, K=200
    )
    mm = materials.muesliMaterial(
        materials.MuesliFiniteStrain.yeohMaterial,
        C=[100, 200, 300],
        K=400,
        compressible=False,
    )

    mm.printDescription()


if __name__ == "__main__":
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

    checkStrainTransformation()
    checkVoigtTransformations()
    checkMaterialConstructors()

    instantiateMuesliMaterials()
