import ikarus as iks
from ikarus import materials

import numpy as np

Emod = 1000
nu = 0.0

svk = materials.StVenantKirchhoff(E=Emod, nu=nu)
ps = svk.asPlaneStress()

E_ = np.array([[-0.375,   0,    0], [0, 0, 0], [0, 0, 0]], dtype=float)
 

print(type(E_))

E = materials.toVoigt(E_)

C_ = materials.tramsformStrain(
    materials.StrainTags.greenLagrangian,
    materials.StrainTags.rightCauchyGreenTensor,
    E_,
)
C = materials.toVoigt(C_)
print(C)

nh = materials.NeoHooke(E=Emod, nu=nu)
psNH = nh.asPlaneStress()
planeStrainNH = nh.asPlaneStrain()


S = svk.stresses(materials.StrainTags.greenLagrangian, E)
print(materials.fromVoigt(S, False))

S = ps.stresses(materials.StrainTags.greenLagrangian, E)
print(materials.fromVoigt(S, False))

S = nh.stresses(materials.StrainTags.rightCauchyGreenTensor, C)
print(materials.fromVoigt(S, False))

S = psNH.stresses(materials.StrainTags.rightCauchyGreenTensor, C)
print(materials.fromVoigt(S, False))

S = planeStrainNH.stresses(materials.StrainTags.rightCauchyGreenTensor, C)
print(materials.fromVoigt(S, False))

# C = svk.tangentModuli(materials.StrainTags.greenLagrangian, E)
# print(C)
