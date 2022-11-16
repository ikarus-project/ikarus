/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#include "physicsHelper.hh"
namespace Ikarus {
  ConvertLameConstants<YoungsModulusAndPoissonsRatio> convertLameConstants(const YoungsModulusAndPoissonsRatio& p_vp) {
    return {p_vp};
  }
  ConvertLameConstants<YoungsModulusAndShearModulus> convertLameConstants(const YoungsModulusAndShearModulus& p_vp) {
    return {p_vp};
  }

  ConvertLameConstants<YoungsModulusAndBulkModulus> convertLameConstants(const YoungsModulusAndBulkModulus& p_vp) {
    return {p_vp};
  }

  ConvertLameConstants<BulkModulusAndLamesFirstParameter> convertLameConstants(
      const BulkModulusAndLamesFirstParameter& p_vp) {
    return {p_vp};
  }

  Eigen::Matrix3d planeStressLinearElasticMaterialTangent(double E, double nu) {
    Eigen::Matrix3d C;
    C.setZero();
    C(0, 0) = C(1, 1) = 1;
    C(0, 1) = C(1, 0) = nu;
    C(2, 2)           = (1 - nu) / 2;
    C *= E / (1 - nu * nu);
    return C;
  }

  Eigen::Matrix<double, 6, 6> linearElasticMaterialTangent3D(double E, double nu) {
    Eigen::Matrix<double, 6, 6> C;
    C.setZero();
    C(0, 0) = C(1, 1) = C(2, 2) = 1 - nu;
    C(0, 1) = C(1, 0) = C(2, 0) = C(0, 2) = C(1, 2) = C(2, 1) = nu;
    C(3, 3) = C(4, 4) = C(5, 5) = (1 - 2 * nu) / 2;
    C *= E / ((1 + nu) * (1 - 2 * nu));
    return C;
  }
}  // namespace Ikarus