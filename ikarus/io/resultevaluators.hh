// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/math.hh>

#include <ikarus/finiteelements/ferequirements.hh>

namespace Dune {
  // Forward declaration
  template <typename ScalarType, int size>
  class FieldVector;
}  // namespace Dune

namespace Ikarus::ResultEvaluators {

  struct VonMises {
    template <typename ElementType, typename FERequirements, int size, typename ScalarType>
    double operator()(const ElementType &fe, const Ikarus::ResultRequirements<FERequirements> &req,
                      const Dune::FieldVector<ScalarType, size> &pos, [[maybe_unused]] int comp) const
        requires(size == 2) {
      Ikarus::ResultTypeMap<ScalarType> res_;
      fe.calculateAt(req, pos, res_);

      const auto &[resultType, sigma] = res_.getSingleResult();
      assert(resultType == ResultType::cauchyStress or resultType == ResultType::PK2Stress);
      const auto s_x  = sigma(0, 0);
      const auto s_y  = sigma(1, 0);
      const auto s_xy = sigma(2, 0);

      return std::sqrt(std::pow(s_x, 2) + Dune::power(s_y, 2) - s_x * s_y + 3 * Dune::power(s_xy, 2));
    }
    static std::string name() { return "VonMises"; }
    static int ncomps() { return 1; }
  };

  struct PrincipalStress {
    template <typename ElementType, typename FERequirements, int size, typename ScalarType>
    double operator()(const ElementType &fe, const Ikarus::ResultRequirements<FERequirements> &req,
                      const Dune::FieldVector<ScalarType, size> &pos, int comp) const requires(size == 2) {
      Ikarus::ResultTypeMap<ScalarType> res_;
      fe.calculateAt(req, pos, res_);

      const auto &[resultType, sigma] = res_.getSingleResult();
      assert(resultType == ResultType::cauchyStress or resultType == ResultType::PK2Stress);
      const auto s_x  = sigma(0, 0);
      const auto s_y  = sigma(1, 0);
      const auto s_xy = sigma(2, 0);

      auto t1 = (s_x + s_y) / 2;
      auto t2 = std::sqrt(Dune::power((s_x - s_y) / 2, 2) + Dune::power(s_xy, 2));

      return comp == 0 ? t1 + t2 : t1 - t2;
    }
    static std::string name() { return "PrincipalStress"; }
    static int ncomps() { return 2; }
  };
}  // namespace Ikarus::ResultEvaluators
