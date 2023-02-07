// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

namespace Ikarus
{
  // clang-format off
  enum class ScalarAffordances {
    noAffordance,
    mechanicalPotentialEnergy,
    microMagneticPotentialEnergy
  };

  enum class VectorAffordances {
    noAffordance,
    forces,
    microMagneticForces
  };

  enum class MatrixAffordances {
    noAffordance,
    stiffness,
    materialstiffness,
    geometricstiffness,
    stiffnessdiffBucklingVector,
    microMagneticHessian,
    mass
  };

  enum class FEParameter {
    noParameter,
    loadfactor,
    time
  };

  enum class FESolutions {
    noSolution,
    displacement,
    velocity,
    director,
    magnetizationAndVectorPotential
  };

  enum class ResultType {
    noType,
    magnetization,
    gradientNormOfMagnetization,
    vectorPotential,
    divergenceOfVectorPotential,
    BField,
    HField,
    cauchyStress,
    director
  };
  // clang-format on
}
