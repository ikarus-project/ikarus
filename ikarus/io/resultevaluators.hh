// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file resultevaluators.hh
 * \brief Ikarus Result Evaluators for special stress quantities
 * \ingroup resultevaluators
 *
 */

#pragma once

#include <dune/common/math.hh>

#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::ResultEvaluators {

/**
 * \brief Struct for calculating von Mises stress
 * \ingroup resultevaluators
 * \details The VonMises struct provides a function call operator to calculate von Mises stress.
 * In 2D, this assumes a plane stress state
 */
struct VonMises
{
  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenMatrix containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return von Mises stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const int comp) const {
    const auto s_x = resultArray(0, 0);
    const auto s_y = resultArray(1, 0);
    if constexpr (R::CompileTimeTraits::RowsAtCompileTime == 3) {
      const auto s_xy = resultArray(2, 0);
      return std::sqrt(Dune::power(s_x, 2) + Dune::power(s_y, 2) - s_x * s_y + 3 * Dune::power(s_xy, 2));
    } else {
      const auto s_z  = resultArray(2, 0);
      const auto s_yz = resultArray(3, 0);
      const auto s_xz = resultArray(4, 0);
      const auto s_xy = resultArray(5, 0);

      return std::sqrt(Dune::power(s_x, 2) + Dune::power(s_y, 2) + Dune::power(s_z, 2) - s_x * s_y - s_x * s_z -
                       s_y * s_z + 3 * (Dune::power(s_xy, 2) + Dune::power(s_xz, 2) + Dune::power(s_yz, 2)));
    }
  }

  /**
   * \brief Get the name of the result type (VonMises)
   * \return String representing the name
   */
  static std::string name() { return "VonMises"; }

  /**
   * \brief Get the number of components in the result (always 1 for VonMises)
   * \return Number of components
   */
  static int ncomps() { return 1; }
};

/**
 * \brief Struct for calculating von Mises stress
 * \ingroup resultevaluators
 * \details The VonMises struct provides a function call operator to calculate von Mises stress.
 * In 2D, this assumes a plane stress state
 */
struct HydrostaticStress
{
  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenMatrix containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return Hydrostatic stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const int comp) const {
    static constexpr int dim = R::CompileTimeTraits::RowsAtCompileTime;
    const auto sigma         = fromVoigt(resultArray, false);
    return 1.0 / dim * sigma.trace();
  }

  /**
   * \brief Get the name of the result type (VonMises)
   * \return String representing the name
   */
  static std::string name() { return "HydrostaticStress"; }

  /**
   * \brief Get the number of components in the result (always 1 for VonMises)
   * \return Number of components
   */
  static int ncomps() { return 1; }
};

/**
 * \brief Struct for calculating principal stresses
 * \ingroup resultevaluators
 * \details The PrincipalStress struct provides a function call operator to calculate principal stresses.
 * The components are ordered in a descending manner ($\sigma_1 > \sigma_2$)
 * \tparam dim dimension of stress state
 */
template <int dim>
requires(dim == 2 or dim == 3)
struct PrincipalStress
{
  /**
   * \brief Calculate the result quantity (principal stress)
   * \param resultArray EigenMatrix containing the stress state in Voigt notation
   * \param comp component of result
   * \return principal stress
   */
  double operator()(const auto& resultArray, const int comp) const {
    auto mat = fromVoigt(resultArray, false);
    Eigen::SelfAdjointEigenSolver<decltype(mat)> eigensolver(mat, Eigen::EigenvaluesOnly);
    return eigensolver.eigenvalues()[dim - 1 - comp];
  }

  /**
   * \brief Get the name of the result type (PrincipalStress)
   * \return String representing the name
   */
  static std::string name() { return "PrincipalStress"; }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  static int ncomps() { return dim; }
};

/**
 * \brief Struct for calculating Triaxiality stresses
 * \ingroup resultevaluators
 * \details The VonMises struct provides a function call operator to calculate von Mises stress.
 * In 2D, this assumes a plane stress state
 */
struct Triaxiality
{
  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenMatrix containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return Triaxiality stress
   */
  template <typename R>
  double operator()(const R& resultArray, const int comp) const {
    auto sigeq = VonMises{}(resultArray, 0);
    auto sigm  = HydrostaticStress{}(resultArray, 0);
    return sigm / sigeq;
  }
  /**
   * \brief Get the name of the result type (PrincipalStress)
   * \return String representing the name
   */
  static std::string name() { return "Triaxiality"; }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  static int ncomps() { return 1; }
};

template <typename ResultEvaluator>
struct PlaneStrainWrapper
{
  template <typename RE>
  PlaneStrainWrapper(RE&& resultEvaluator, double nu)
      : underlying_(std::forward<RE>(resultEvaluator)),
        nu_(nu) {}

  template <typename R>
  double operator()(const R& resultArray, const int comp) const {
    static_assert(R::CompileTimeTraits::RowsAtCompileTime == 3, "PlaneStrainWrapper is only valid for 2D.");
    auto sigZ                     = nu_ * (resultArray[0] + resultArray[1]);
    auto enlargedResultArray      = Eigen::Vector<double, 6>::Zero().eval();
    enlargedResultArray.head<2>() = resultArray.template head<2>();
    enlargedResultArray[3]        = sigZ;
    enlargedResultArray[5]        = resultArray[2];

    return underlying_(enlargedResultArray, comp);
  }
  /**
   * \brief Get the name of the result type (PrincipalStress)
   * \return String representing the name
   */
  static std::string name() { return ResultEvaluator::name(); }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  static int ncomps() { ResultEvaluator::ncomps(); }

private:
  ResultEvaluator underlying_;
  double nu_;
};

template <typename ResultEvaluator>
PlaneStrainWrapper(ResultEvaluator&&, double) -> PlaneStrainWrapper<ResultEvaluator>;

} // namespace Ikarus::ResultEvaluators
