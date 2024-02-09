// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later


#pragma once

#include <ikarus/utils/math.hh>
/**
*
* \ingroup FEParameterTags
* \brief A strongly typed enum class representing the type of the result request
*/


namespace Ikarus::ResultType
{
#define REGISTER_RT(structName) \
friend auto toString(structName){return #structName;}

namespace Impl {
  template <int dim>
  constexpr int matrixSize() {
    return (-1 + Ikarus::ct_sqrt(1 + 8 * dim)) / 2;
  }
}

struct noType;
struct magnetization;
struct gradientNormOfMagnetization;
struct vectorPotential;
struct divergenceOfVectorPotential;
struct BField;
struct HField;
struct cauchyStress;
struct PK2Stress;
struct linearStress;
struct director;

struct noType
{
  REGISTER_RT(noType);
};

struct magnetization
{
  REGISTER_RT(magnetization);
};

struct gradientNormOfMagnetization
{
  REGISTER_RT(gradientNormOfMagnetization);
};
struct vectorPotential
{
  REGISTER_RT(vectorPotential);
};
struct divergenceOfVectorPotential
{
  REGISTER_RT(divergenceOfVectorPotential);
};
struct BField
{
  REGISTER_RT(BField);
};
struct HField
{
  REGISTER_RT(HField);
};
struct cauchyStress
{
  REGISTER_RT(cauchyStress);
};
struct PK2Stress
{
  REGISTER_RT(PK2Stress);
};
struct linearStress
{
  REGISTER_RT(linearStress);

};
struct director
{
  REGISTER_RT(director);

};

}