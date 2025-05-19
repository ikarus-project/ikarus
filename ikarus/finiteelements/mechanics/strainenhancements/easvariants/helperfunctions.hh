// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file helperfunctions.hh
 * \brief Definition of the certain helper functions used to transform strain measures.
 *
 * \ingroup eas
 */

#pragma once

#include <dune/localfefunctions/eigenDuneTransformations.hh>

#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS::Impl {
template <typename GEO>
auto transformationMatrixAtCenterWithDetJ(const GEO& geometry) {
  const auto& referenceElement = Dune::ReferenceElements<double, GEO::mydimension>::general(geometry.type());
  const auto quadPos0          = referenceElement.position(0, 0);
  const auto detJ0             = geometry.integrationElement(quadPos0);

  return (transformationMatrix(geometry, quadPos0) * detJ0).eval();
}

template <typename GEO>
void transformDisplacementGradient(const GEO& geometry, Eigen::Matrix<double, GEO::mydimension, GEO::mydimension>& H,
                                   const Dune::FieldVector<double, GEO::mydimension>& quadPos) {
  const auto& referenceElement = Dune::ReferenceElements<double, GEO::mydimension>::general(geometry.type());
  const auto quadPos0          = referenceElement.position(0, 0);
  const auto detJ0             = geometry.integrationElement(quadPos0);
  const auto detJ              = geometry.integrationElement(quadPos);
  const auto JInv              = Dune::toEigen(geometry.jacobianInverseTransposed(quadPos0)).eval(); // pre-transposed
  H                            = ((detJ0 / detJ) * JInv * H * JInv.transpose()).eval();
}
} // namespace Ikarus::EAS::Impl