//
// Created by lex on 07/03/2022.
//

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
}  // namespace Ikarus