# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

from dune.common.hashit import hashIt
from ikarus.generator import MySimpleGenerator

def volumeLoad(f,d):
    generator = MySimpleGenerator("VolumeLoadPre", "Ikarus::Python")
    element_type = f"Ikarus::VolumeLoadPre<{d}>"

    includes = []
    includes += ["ikarus/python/finiteelements/loads/volume.hh"]
    moduleName = "VolumeLoadPre_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.VolumeLoadPre( f)

def volumeLoad2D(f):
    return volumeLoad(f,2)

def volumeLoad3D(f):
    return volumeLoad(f,3)

def neumannBoundaryLoad(boundaryPatch, f):
    generator = MySimpleGenerator("NeumannBoundaryLoadPre", "Ikarus::Python")
    element_type = f"Ikarus::NeumannBoundaryLoadPre<{boundaryPatch.gridView().cppTypeName}>"

    includes = []
    includes += ["ikarus/python/finiteelements/loads/traction.hh"]
    includes += boundaryPatch._includes
    moduleName = "NeumannBoundaryLoadPre_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    return module.NeumannBoundaryLoadPre(boundaryPatch, f)

def registerPreElement(name,includes,element_type,*args):
    generator = MySimpleGenerator(name, "Ikarus::Python")
    includes+= ["ikarus/python/finiteelements/registerpreelement.hh"]

    moduleName = name + "_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    f = getattr(module, name)
    return f(*args)

def nonLinearElastic(mat):
    includes = ["ikarus/finiteelements/mechanics/nonlinearelastic.hh"]

    element_type = f"Ikarus::NonLinearElasticPre<{mat.cppTypeName}>"
    includes += mat._includes
    return  registerPreElement("NonLinearElasticPre",includes,element_type,mat)


def linearElastic(youngs_modulus,nu):
    includes = ["ikarus/finiteelements/mechanics/linearelastic.hh"]
    element_type = "Ikarus::LinearElasticPre"
    return  registerPreElement("LinearElasticPre",includes,element_type,youngs_modulus,nu)

def eas(numberofparameters):
    includes = ["ikarus/finiteelements/mechanics/enhancedassumedstrains.hh"]
    element_type = "Ikarus::EnhancedAssumedStrainsPre"
    return  registerPreElement("EnhancedAssumedStrainsPre",includes,element_type,numberofparameters)

def kirchhoffLoveShell(youngs_modulus,nu, thickness):
    includes = ["ikarus/finiteelements/mechanics/kirchhoffloveshell.hh"]
    element_type = "Ikarus::KirchhoffLoveShellPre"
    return  registerPreElement("KirchhoffLoveShellPre",includes,element_type,youngs_modulus,nu,thickness)



def makeFE(basis,*skills):

    includes = ["ikarus/python/finiteelements/fe.hh"]
    preFE_type = f"Ikarus::PreFE<{basis.cppTypeName},Ikarus::FERequirements<Eigen::Ref<Eigen::VectorXd>>,true,true>"
    element_type = f"Ikarus::FE<{preFE_type},"

    for arg in skills:
        element_type += arg.cppTypeName+"::Skill,"
    element_type = element_type[:-1]
    element_type+=">"

    generator = MySimpleGenerator("FE", "Ikarus::Python")

    includes += basis._includes
    for arg in skills:
        includes += arg._includes

    moduleName = "FE" + "_" + hashIt(element_type)
    module = generator.load(
        includes=includes, typeName=element_type, moduleName=moduleName
    )

    tupleskills = tuple(skills)
    return module.FE(basis,tupleskills)
