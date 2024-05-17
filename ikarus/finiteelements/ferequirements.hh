// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file ferequirements.hh
 * \brief Definition of the LinearElastic class for finite element mechanics computations.
 * \ingroup finiteelements
 */

#pragma once

#include <iosfwd>
#include <map>
#include <set>
#include <vector>

#include <dune/common/exceptions.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus {
// clang-format off

  /**
 * ScalarAffordance
 * \ingroup Affordancetags
 * \brief A strongly typed enum class representing the scalar affordance
 */
// cppcheck-suppress MAKE_ENUM
 MAKE_ENUM(ScalarAffordance,
            noAffordance,
            mechanicalPotentialEnergy,
            microMagneticPotentialEnergy
      );

  /**
* VectorAffordance
* \ingroup Affordancetags
* \brief A strongly typed enum class representing the vector affordance
*/
  MAKE_ENUM(VectorAffordance,
            noAffordance,
            forces,
            microMagneticForces
      );

  /**
* MatrixAffordance
* \ingroup Affordancetags
* \brief A strongly typed enum class representing the matrix affordance
*/
  MAKE_ENUM(MatrixAffordance,
            noAffordance,
            stiffness,
            materialstiffness,
            geometricstiffness,
            stiffnessdiffBucklingVector,
            microMagneticHessian,
            mass
      );

  /**
*
* \ingroup FEParameterTags
* \brief A strongly typed enum class representing the FE parameter
*/
  MAKE_ENUM(FEParameter,
            noParameter,
            loadfactor,
            time
      );

  /**
*
* \ingroup FEParameterTags
* \brief A strongly typed enum class representing the type of the solutions vectors
*/
  MAKE_ENUM(FESolutions,
            noSolution,
            displacement,
            velocity,
            director,
            magnetizationAndVectorPotential
      );

// clang-format on

/**
 * \brief Concept to check if a given type is one of the predefined affordance enums or the AffordanceCollection
 */
template <typename T>
concept FEAffordance = std::is_same_v<std::remove_cvref_t<T>, ScalarAffordance> or
                       std::is_same_v<std::remove_cvref_t<T>, VectorAffordance> or
                       std::is_same_v<std::remove_cvref_t<T>, MatrixAffordance>;

/**
 * \brief Struct representing a collection of affordances.
 */
template <FEAffordance... Affos>
requires(sizeof...(Affos) <= 3)
struct AffordanceCollection : public std::tuple<Affos...>
{
  using Base = std::tuple<Affos...>;

  AffordanceCollection() = default;
  constexpr AffordanceCollection(Affos... affos)
      : Base(affos...) {}

  static constexpr bool hasScalarAffordance = traits::hasType<ScalarAffordance, std::tuple<Affos...>>::value;
  static constexpr bool hasVectorAffordance = traits::hasType<VectorAffordance, std::tuple<Affos...>>::value;
  static constexpr bool hasMatrixAffordance = traits::hasType<MatrixAffordance, std::tuple<Affos...>>::value;
  /**
   * \brief Check if a specific affordance is present in the requirements.
   *
   * This function checks if the specified affordance is present in the requirements.
   *
   * \tparam Affordance Type of affordance to be checked.
   * \param affordance The affordance to be checked.
   * \return True if the affordance is present, false otherwise.
   */
  template <FEAffordance Affordance>
  constexpr bool hasAffordance(Affordance&& affordances) const {
    if constexpr (std::is_same_v<std::remove_cvref_t<Affordance>, AffordanceCollection>)
      return affordances == *this;
    else {
      if constexpr (std::is_same_v<std::remove_cvref_t<Affordance>, ScalarAffordance>) {
        if constexpr (hasScalarAffordance)
          return affordances == std::get<ScalarAffordance>(*this);
        else
          return false;
      }
      if constexpr (std::is_same_v<std::remove_cvref_t<Affordance>, VectorAffordance>) {
        if constexpr (hasVectorAffordance)
          return affordances == std::get<VectorAffordance>(*this);
        else
          return false;
      }
      if constexpr (std::is_same_v<std::remove_cvref_t<Affordance>, MatrixAffordance>) {
        if constexpr (hasMatrixAffordance)
          return affordances == std::get<MatrixAffordance>(*this);
        else
          return false;
      }
    }
    return false;
  }

  auto scalarAffordance() const
  requires hasScalarAffordance
  {
    return std::get<ScalarAffordance>(*this);
  }
  auto vectorAffordance() const
  requires hasVectorAffordance
  {
    return std::get<VectorAffordance>(*this);
  }
  auto matrixAffordance() const
  requires hasMatrixAffordance
  {
    return std::get<MatrixAffordance>(*this);
  }
};

inline constexpr VectorAffordance forces = VectorAffordance::forces;

inline constexpr MatrixAffordance stiffness                   = MatrixAffordance::stiffness;
inline constexpr MatrixAffordance stiffnessdiffBucklingVector = MatrixAffordance::stiffnessdiffBucklingVector;
inline constexpr MatrixAffordance mass                        = MatrixAffordance::mass;
inline constexpr ScalarAffordance potentialEnergy             = ScalarAffordance::mechanicalPotentialEnergy;

auto vectorAffordance(MatrixAffordance affoM) {
  if (affoM == MatrixAffordance::stiffness)
    return VectorAffordance::forces;
  else if (affoM == MatrixAffordance::microMagneticHessian)
    return VectorAffordance::microMagneticForces;
  else
    return VectorAffordance::noAffordance;
}

auto scalarAffordance(MatrixAffordance affoM) {
  if (affoM == MatrixAffordance::stiffness)
    return ScalarAffordance::mechanicalPotentialEnergy;
  else if (affoM == MatrixAffordance::microMagneticHessian)
    return ScalarAffordance::microMagneticPotentialEnergy;
  else
    return ScalarAffordance::noAffordance;
}

auto scalarAffordance(VectorAffordance affoV) {
  if (affoV == VectorAffordance::forces)
    return ScalarAffordance::mechanicalPotentialEnergy;
  else if (affoV == VectorAffordance::microMagneticForces)
    return ScalarAffordance::microMagneticPotentialEnergy;
  else
    return ScalarAffordance::noAffordance;
}

namespace AffordanceCollections {
  inline constexpr AffordanceCollection elastoStatics(ScalarAffordance::mechanicalPotentialEnergy,
                                                      VectorAffordance::forces, MatrixAffordance::stiffness);
}

namespace Impl {
  template <typename T>
  struct DeduceRawVectorType
  {
    using Type = T;
  };

  template <typename T>
  struct DeduceRawVectorType<std::reference_wrapper<T>>
  {
    using Type = T;
  };

  template <typename T>
  struct DeduceRawVectorType<Eigen::Ref<T>>
  {
    using Type = Eigen::Ref<T>;
  };
} // namespace Impl

/**
 * \class FERequirements
 * \brief Class representing the requirements for finite element calculations.
 *
 * This class defines the requirements for finite element calculations, including the types of solution vectors
 * and parameters needed. It provides methods to add affordances, insert parameters, and manage global solution
 * vectors.
 *
 * \tparam SV Type of the solution vector, defaulting to std::reference_wrapper<Eigen::VectorXd>.
 * \tparam PM Type of the parameter, defaulting to std::reference_wrapper<double>.
 *
 */
template <FESolutions sol, FEParameter para, typename SV = Eigen::VectorXd, typename PM = double>
class FERequirements
{
public:
  static constexpr FESolutions globalSolutionTag = sol;
  static constexpr FEParameter parameterTag      = para;
  using SolutionVectorType                       = SV;
  // using SolutionVectorTypeRaw = typename Impl::DeduceRawVectorType<std::remove_cvref_t<SV>>::Type;
  using ParameterType = PM;

  FERequirements() = default;
  FERequirements(SolutionVectorType& solVec, ParameterType& parameter)
      : sol_{solVec},
        parameter_{parameter} {}

  /**
   * \brief Insert a parameter into the requirements.
   *
   * This function inserts the specified parameter into the requirements.
   *
   * \param val Reference to the raw parameter value.
   * \return Reference to the updated FERequirements instance.
   */
  FERequirements& insertParameter(ParameterType& val) {
    parameter_ = val;
    return *this;
  }

  /**
   * \brief Insert a global solution vector into the requirements.
   *
   * This function inserts the specified global solution vector into the requirements.
   *
   * \param solVec Reference to the raw global solution vector.
   * \return Reference to the updated FERequirements instance.
   */
  FERequirements& insertGlobalSolution(SolutionVectorType& solVec) {
    sol_ = solVec;
    return *this;
  }

  /**
   * \brief Get the global solution vector.
   *
   * \return Reference to the raw global solution vector.
   *
   */
  SolutionVectorType& globalSolution() { return sol_.value().get(); }

  /**
   * \brief Get the  global solution vector.
   *   *
   * \return Const reference to the global solution vector.
   *
   */
  const SolutionVectorType& globalSolution() const { return sol_.value().get(); }

  /**
   * \brief Get the parameter value.
   *
   *
   * \return Const reference to the parameter value.
   *
   */
  const ParameterType& parameter() const {
    return parameter_.value().get();
    ;
  }

  /**
   * \brief Get the parameter value.
   *
   * \return Reference to the parameter value.
   *
   */
  ParameterType& parameter() { return parameter_.value().get(); }

  /**
   * \brief Tells if the class contains all needed values
   *
   * \return Bool if the all values are assigned
   *
   */
  bool populated() const { return sol_.has_value() and parameter_.has_value(); }

private:
  std::optional<std::reference_wrapper<SolutionVectorType>> sol_;
  std::optional<std::reference_wrapper<ParameterType>> parameter_;
};

template <FESolutions sol, FEParameter para, bool wrapWithRef = false, typename SV = Eigen::VectorXd,
          typename PM = double>
struct FERequirementsFactory
{
private:
  using typeEigen = std::conditional_t<wrapWithRef and Ikarus::Concepts::EigenMatrix<SV>, Eigen::Ref<SV>, SV>;

public:
  using type = FERequirements<sol, para, typeEigen, PM>;
};

} // namespace Ikarus
