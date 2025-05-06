// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
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
#include <ikarus/utils/basis.hh>
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
template <FEAffordance... Affordances>
requires(sizeof...(Affordances) <= 3)
struct AffordanceCollection : public std::tuple<Affordances...>
{
  using Base = std::tuple<Affordances...>;

  AffordanceCollection() = default;
  constexpr AffordanceCollection(Affordances... affordances)
      : Base(affordances...) {}

  static constexpr bool hasScalarAffordance = traits::hasType<ScalarAffordance, std::tuple<Affordances...>>::value;
  static constexpr bool hasVectorAffordance = traits::hasType<VectorAffordance, std::tuple<Affordances...>>::value;
  static constexpr bool hasMatrixAffordance = traits::hasType<MatrixAffordance, std::tuple<Affordances...>>::value;
  /**
   * \brief Check if a specific affordance is present in the requirements.
   *
   * This function checks if the specified affordance is present in the requirements.
   *
   * \tparam Affordance Type of affordance to be checked.
   * \param affordance The affordance to be checked. This can also be an affordance collection.
   * \return True if the affordance is present, false otherwise.
   */
  template <FEAffordance Affordance>
  constexpr bool hasAffordance(Affordance&& affordances) const {
    using AffordanceRaw = std::remove_cvref_t<Affordance>;
    if constexpr (std::is_same_v<AffordanceRaw, AffordanceCollection>)
      return affordances == *this;
    else {
      if constexpr (std::is_same_v<AffordanceRaw, ScalarAffordance>) {
        if constexpr (hasScalarAffordance)
          return affordances == std::get<ScalarAffordance>(*this);
        else
          return false;
      }
      if constexpr (std::is_same_v<AffordanceRaw, VectorAffordance>) {
        if constexpr (hasVectorAffordance)
          return affordances == std::get<VectorAffordance>(*this);
        else
          return false;
      }
      if constexpr (std::is_same_v<AffordanceRaw, MatrixAffordance>) {
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

auto vectorAffordance(MatrixAffordance affordanceM) {
  if (affordanceM == MatrixAffordance::stiffness)
    return VectorAffordance::forces;
  else if (affordanceM == MatrixAffordance::microMagneticHessian)
    return VectorAffordance::microMagneticForces;
  else
    return VectorAffordance::noAffordance;
}

auto scalarAffordance(MatrixAffordance affordanceM) {
  if (affordanceM == MatrixAffordance::stiffness)
    return ScalarAffordance::mechanicalPotentialEnergy;
  else if (affordanceM == MatrixAffordance::microMagneticHessian)
    return ScalarAffordance::microMagneticPotentialEnergy;
  else
    return ScalarAffordance::noAffordance;
}

auto scalarAffordance(VectorAffordance affordanceV) {
  if (affordanceV == VectorAffordance::forces)
    return ScalarAffordance::mechanicalPotentialEnergy;
  else if (affordanceV == VectorAffordance::microMagneticForces)
    return ScalarAffordance::microMagneticPotentialEnergy;
  else
    return ScalarAffordance::noAffordance;
}

namespace AffordanceCollections {
  inline constexpr AffordanceCollection elastoStatics(ScalarAffordance::mechanicalPotentialEnergy,
                                                      VectorAffordance::forces, MatrixAffordance::stiffness);
} // namespace AffordanceCollections

/**
 * \class FERequirements
 * \brief Class representing the requirements for finite element calculations.
 *
 * This class defines the requirements for finite element calculations, including the types of solution vectors
 * and parameters needed.
 *
 * \tparam sol The finite element solution tag.
 * \tparam para The finite element parameter tag.
 * \tparam SV Type of the solution vector, defaulting to Eigen::VectorXd.
 * \tparam PM Type of the parameter, defaulting to double.
 */
template <FESolutions sol, FEParameter para, typename SV = Eigen::VectorXd, typename PM = double>
class FERequirements
{
public:
  using SolutionVectorType = SV; ///< Type of the solution vector.
  using ParameterType      = PM; ///< Type of the parameter.

  /**
   * \brief Default constructor.
   */
  FERequirements() = default;

  /**
   * \brief Constructor initializing the solution vector and parameter.
   *
   * \param solVec Solution vector.
   * \param parameter Parameter value.
   */
  template <typename SV2 = SV, typename PM2 = PM>
  FERequirements(SV2&& solVec, PM2&& parameter)
      : sol_(std::make_unique<SV>(std::forward<SV2>(solVec))),
        parameter_(std::make_unique<PM>(std::forward<PM2>(parameter))) {}

  /**
   * \brief Constructor from a basis
   * \tparam PB The type of the basis.
   * \param basis the basis
   */
  template <typename PB>
  FERequirements(const Ikarus::BasisHandler<PB>& basis)
      : sol_(std::make_unique<SV>(SV::Zero(basis.flat().size()))),
        parameter_(std::make_unique<PM>(0.0)) {}

  /**
   * \brief Copy constructor.
   *
   * Performs a deep copy of the solution vector and parameter.
   *
   * \param other The object to copy from.
   */
  FERequirements(const FERequirements& other)
      : sol_(other.sol_ ? std::make_unique<SV>(*other.sol_) : nullptr),
        parameter_(other.parameter_ ? std::make_unique<PM>(*other.parameter_) : nullptr) {}

  /**
   * \brief Copy assignment operator.
   *
   * Performs a deep copy of the solution vector and parameter.
   *
   * \param other The object to assign from.
   * \return Reference to the updated FERequirements instance.
   */
  FERequirements& operator=(const FERequirements& other) {
    if (this != &other) {
      sol_       = other.sol_ ? std::make_unique<SV>(*other.sol_) : nullptr;
      parameter_ = other.parameter_ ? std::make_unique<PM>(*other.parameter_) : nullptr;
    }
    return *this;
  }

  /**
   * \brief Insert a parameter into the requirements.
   *
   * This function inserts the specified parameter into the requirements.
   *
   * \param val Reference to the raw parameter value.
   * \return Reference to the updated FERequirements instance.
   */
  FERequirements& insertParameter(const PM& val) {
    parameter_ = std::make_unique<PM>(val);
    return *this;
  }

  /**
   * \brief Insert a global solution vector into the requirements.
   *
   * This function inserts the specified global solution vector into the requirements and deletes the old one.
   *
   * \param solVec Reference to the raw global solution vector.
   * \return Reference to the updated FERequirements instance.
   */
  template <typename SV2 = SolutionVectorType>
  FERequirements& insertGlobalSolution(SV2&& solVec) {
    sol_ = std::make_unique<SV>(std::forward<SV2>(solVec));
    return *this;
  }

  /**
   * \brief Get the global solution vector.
   *
   * \return Reference to the raw global solution vector.
   */
  const SolutionVectorType& globalSolution() const {
    if (!sol_)
      DUNE_THROW(Dune::InvalidStateException, "Solution vector is not initialized.");
    return *sol_;
  }

  /**
   * \brief Get the global solution vector.
   *
   * \return Reference to the raw global solution vector.
   */
  SV& globalSolution() {
    if (!sol_)
      DUNE_THROW(Dune::InvalidStateException, "Solution vector is not initialized.");

    return *sol_;
  }

  /**
   * \brief Get the parameter value.
   *
   * \return Reference to the parameter value.
   */
  const PM& parameter() const {
    if (!parameter_)
      DUNE_THROW(Dune::InvalidStateException, "Parameter is not initialized.");

    return *parameter_;
  }

  /**
   * \brief Get the parameter value.
   *
   * \return Reference to the parameter value.
   */
  PM& parameter() {
    if (!parameter_)
      DUNE_THROW(Dune::InvalidStateException, "Parameter is not initialized.");

    return *parameter_;
  }

  /**
   * \brief Tells if the class contains all needed values.
   *
   * \return Bool indicating if all values are assigned.
   */
  bool populated() const { return sol_ and parameter_; }

  /**
   * \brief Enables the usage of the class as a solution vector.
   *
   * \tparam T The type of the value to add to the solution vector.
   * \param rhs The value to add to the solution vector.
   * \return Reference to the updated solution vector.
   */
  template <typename T>
  SV& operator+=(const T& rhs) {
    if (!sol_)
      DUNE_THROW(Dune::InvalidStateException, "Solution vector is not initialized.");
    *sol_ += rhs;
    return *sol_;
  }

private:
  std::unique_ptr<SV> sol_;       ///< Unique pointer to the solution vector.
  std::unique_ptr<PM> parameter_; ///< Unique pointer to the parameter.
};

} // namespace Ikarus
