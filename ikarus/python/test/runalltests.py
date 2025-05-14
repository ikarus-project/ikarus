# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import dirichletvaluetest
import linearelastictest
import nonlinearelastictest
import kltest
import testmaterials
import trusstest
import vtkwritertest

print("Running dirichletvaluetest")
dirichletvaluetest.main()

print("Running linearelastictest")
linearelastictest.main()

print("Running nonlinearelastictest")
nonlinearelastictest.main()

print("Running testmaterials")
testmaterials.main()

print("Running trusstest")
trusstest.main()

print("Running vtkwritertest")
vtkwritertest.main()

# Currently disabled because of brakage of dune-iga
# kltest.main()
