// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file fefactory.hh
 * \brief Contains the definition of the FEFactory class
 * \ingroup finiteelements
 */

#pragma once
#include <ikarus/finiteelements/febase.hh>
namespace Ikarus {
/**
 * \brief FEFactory is a convenient wrapper to forward arguments to PreFE and create and
 * construct a factory of finite elements.
 * \tparam BH The type of the basis handler.
 * \tparam SK Type of the skills.
 * \tparam useFlat A boolean indicating if the type of the underlying basis is of the flat or the untouched version.
 * \tparam FER The requirements for the finite element.
 * \tparam useEigenRef A boolean flag indicating whether to use Eigen references.
 */
template <typename BH, typename SK, typename FER = FERequirements<>, bool useFlat = true, bool useEigenRef = false>
struct FEFactory
{
  using Skills = SK;

private:
  const BH* basisHandler_;
  SK skills;

public:
  /**
   * \brief Move constructor for FEFactory
   * \param basisHandler The basis handler.
   * \param sk Skill arguments.
   */
  FEFactory(const BH& basisHandler, SK&& sk)
      : basisHandler_{&basisHandler},
        skills{std::move(sk)} {}
  /**
   * \brief Constructor for FEFactory
   * \param basisHandler The basis handler.
   * \param sk Skill arguments.
   */
  FEFactory(const BH& basisHandler, const SK& sk)
      : basisHandler_{&basisHandler},
        skills{sk} {}

  auto operator()() {
    return std::apply(
        [&]<typename... Args>(Args&&... args) {
          // the template would not be needed,
          // when https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p1814r0.html
          // will be implemented in clang. It is already implemented in gcc 12.2
          typename PreFE<BH, FER, useFlat, useEigenRef>::template FE<std::decay_t<Args>::template Skill...> fe(
              *basisHandler_, std::forward<Args>(args)...);

          return fe;
        },
        skills.args);
  }
};

/**
 * \brief A function to create a finite element using the flat version of the basis.
 * \tparam BH The type of the basis handler.
 * \tparam SK Type of the skills.
 * \tparam useFlat A boolean indicating if the  the underlying basis should be handed out as flat or the untouched
 * version. \tparam FER The requirements for the finite element. \tparam useEigenRef A boolean flag indicating whether
 * to use Eigen references. \param basisHandler The basis handler. \param sk Skill arguments. \return An FEFactory
 * object.
 */
template <typename FER = FERequirements<>, bool useFlat = true, bool useEigenRef = false, typename BH, typename SK>
auto makeFE(const BH& basisHandler, SK&& sk) {
  FEFactory<BH, SK, FER, useFlat, useEigenRef> factory(basisHandler, std::forward<SK>(sk));

  return factory();
}

/**
 * \brief A function to create a finite element using the untouched version of the basis.
 * \tparam BH The type of the basis handler.
 * \tparam SK Type of the skills.
 * \tparam FER The requirements for the finite element.
 * \tparam useEigenRef A boolean flag indicating whether to use Eigen references.
 * \param basisHandler The basis handler.
 * \param sk Skill arguments.
 * \return An FEFactory object.
 */
template <typename FER = FERequirements<>, bool useEigenRef = false, typename BH, typename SK>
auto makeFEWithUnTouchedBasis(const BH& basisHandler, SK&& sk) {
  FEFactory<BH, SK, FER, false, useEigenRef> factory(basisHandler, std::forward<SK>(sk));

  return factory();
}

} // namespace Ikarus
