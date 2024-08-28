import ikarus as iks
from ikarus import materials

import numpy as np

svk = materials.StVenantKirchhoff(E=1000, nu=0.3)
ps = svk.asPlaneStress()

E_ = np.array([[100, 0, 1], [0, 230, 110], [1, 110, 50]], dtype=float)

print(type(E_))

E = materials.toVoigt(E_)

C_ = materials.tramsformStrain(
    materials.StrainTags.greenLagrangian,
    materials.StrainTags.rightCauchyGreenTensor,
    E_,
)
C = materials.toVoigt(C_)
print(C)


S = svk.stresses(materials.StrainTags.greenLagrangian, E)
print(materials.fromVoigt(S))

S = ps.stresses(materials.StrainTags.greenLagrangian, E)
print(materials.fromVoigt(S))

nh = materials.NeoHooke(E=1000, nu=0.3)
psNH = nh.asPlaneStress()

S = nh.stresses(materials.StrainTags.greenLagrangian, C)
print(materials.fromVoigt(S))

S = psNH.stresses(materials.StrainTags.greenLagrangian, C)
print(materials.fromVoigt(S))


# C = svk.tangentModuli(materials.StrainTags.greenLagrangian, E)
# print(C)
