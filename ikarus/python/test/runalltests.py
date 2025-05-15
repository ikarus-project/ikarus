# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import runpy

print("Running dirichletvaluetest")
runpy.run_module("dirichletvaluetest", run_name="__main__")

print("Running linearelastictest")
runpy.run_module("linearelastictest", run_name="__main__")

print("Running nonlinearelastictest")
runpy.run_module("nonlinearelastictest", run_name="__main__")

print("Running testmaterials")
runpy.run_module("testmaterials", run_name="__main__")

print("Running trusstest")
runpy.run_module("trusstest", run_name="__main__")

print("Running vtkwritertest")
runpy.run_module("vtkwritertest", run_name="__main__")

print("Running kltest")
runpy.run_module("kltest", run_name="__main__")
