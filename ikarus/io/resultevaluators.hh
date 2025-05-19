// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file resultevaluators.hh
 * \brief Ikarus Result Evaluators for special stress quantities
 * \ingroup resultevaluators
 *
 */

#pragma once

#include <dune/common/math.hh>

#include <Eigen/Geometry>

#include <ikarus/finiteelements/mechanics/materials.hh>
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
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return von Mises stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const auto& pos, [[maybe_unused]] const auto& fe,
                    [[maybe_unused]] const int comp) const {
    auto sigma = fromVoigt(resultArray, false);

    auto I2 = 1.0 / 2.0 * (sigma.squaredNorm() - 1.0 / 3.0 * Dune::power(sigma.trace(), 2));
    return std::sqrt(3.0 * I2);
  }

  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "VonMises"; }

  /**
   * \brief Get the number of components in the result (always 1 for VonMises)
   * \return Number of components
   */
  constexpr static int ncomps() { return 1; }
};

/**
 * \brief Struct for calculating hydrostatic stress
 * \ingroup resultevaluators
 * \details The HydrostaticStress struct provides a function call operator to calculate hydrostatic stress.
 * In 2D, this assumes a plane stress state. Furthermore the stress components are divided by 2 in 2D.
 */
struct HydrostaticStress
{
  /**
   * \brief Calculate the result quantity (hydrostatic stress).
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return Hydrostatic stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const auto& pos, [[maybe_unused]] const auto& fe,
                    [[maybe_unused]] const int comp) const {
    const auto sigma = fromVoigt(resultArray, false);
    return 1.0 / sigma.rows() * sigma.trace();
  }

  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "HydrostaticStress"; }

  /**
   * \brief Get the number of components in the result (always 1 for hydrostatic stress)
   * \return Number of components
   */
  constexpr static int ncomps() { return 1; }
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
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result
   * \return principal stress
   */
  double operator()(const auto& resultArray, [[maybe_unused]] const auto& pos, [[maybe_unused]] const auto& fe,
                    const int comp) const {
    auto mat = fromVoigt(resultArray, false);
    Eigen::SelfAdjointEigenSolver<decltype(mat)> eigensolver(mat, Eigen::EigenvaluesOnly);
    return eigensolver.eigenvalues().reverse()[comp];
  }

  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "PrincipalStress"; }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  constexpr static int ncomps() { return dim; }
};

/**
 * \brief Struct for calculating stress triaxiality
 * \ingroup resultevaluators
 * \details The Triaxiality struct provides a function call operator to calculate stress triaxiality.
 * In 2D, this assumes a plane stress state
 */
struct Triaxiality
{
  /**
   * \brief Calculate the result quantity (stress triaxiality)
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result (not used here)
   * \tparam R Type of the matrix
   * \return Triaxiality stress
   */
  template <typename R>
  double operator()(const R& resultArray, [[maybe_unused]] const auto& pos, [[maybe_unused]] const auto& fe,
                    [[maybe_unused]] const int comp) const {
    auto sigeq = VonMises{}(resultArray, pos, fe, 0);
    auto sigm  = HydrostaticStress{}(resultArray, pos, fe, 0);
    return sigm / sigeq;
  }
  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "Triaxiality"; }

  /**
   * \brief Get the number of components in the result  (always 1 for stress triaxiality)
   * \return Number of components
   */
  constexpr static int ncomps() { return 1; }
};

/**
 * \brief Struct for calculating the 2d polar stress. The center of the coordinate system is to be passed to the
 * evaluator.
 * \ingroup resultevaluators
 */
struct PolarStress
{
  PolarStress(const Dune::FieldVector<double, 2>& origin)
      : origin_(origin) {}

  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenVector containing the stress state in Voigt notation
   * \param comp component of result
   * \tparam R Type of the matrix
   * \return von Mises stress
   */
  template <typename R>
  double operator()(const R& resultArray, const auto& pos, const auto& fe, const int comp) const {
    static_assert(R::CompileTimeTraits::RowsAtCompileTime == 3, "PolarStress is only valid for 2D.");
    if (comp > 2)
      DUNE_THROW(Dune::RangeError, "PolarStress: Comp out of range.");

    // Offset to center the coordinate system in the reference geometry
    Dune::FieldVector<double, 2> posGlobal = fe.geometry().global(pos) - origin_;
    auto theta                             = std::atan2(posGlobal[1], posGlobal[0]);

    const auto sigma = fromVoigt(resultArray, false);
    Eigen::Rotation2D<double> r(theta);

    // deliberately not evaluating this that it stays an expression for below
    auto polarStress = r.inverse() * sigma * r;
    switch (comp) {
      case 0:
        return polarStress(0, 0);
      case 1:
        return polarStress(1, 1);
      case 2:
        return polarStress(1, 0);
      default:
        __builtin_unreachable();
    }
  }

  /**
   * \brief Get the name of the result type
   * \return String representing the name
   */
  constexpr static std::string name() { return "PolarStress"; }

  /**
   * \brief Get the number of components in the result
   * \return Number of components
   */
  constexpr static int ncomps() { return 3; }

private:
  Dune::FieldVector<double, 2> origin_;
};

} // namespace Ikarus::ResultEvaluators
