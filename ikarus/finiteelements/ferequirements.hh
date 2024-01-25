// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file ferequirements.hh
 * @brief Definition of the LinearElastic class for finite element mechanics computations.
 * @ingroup finiteelements
 */

#pragma once

#include <iosfwd>
#include <map>
#include <set>
#include <vector>

#include <dune/common/exceptions.hh>

#include <Eigen/Core>

#include <ikarus/utils/makeenum.hh>

namespace Ikarus {
// clang-format off
  /**
 * ScalarAffordances
 * \ingroup Affordancestags
 * \brief A strongly typed enum class representing the scalar affordance
 */
  MAKE_ENUM(ScalarAffordances,
            noAffordance,
            mechanicalPotentialEnergy,
            microMagneticPotentialEnergy
      );

  /**
* VectorAffordances
* \ingroup Affordancestags
* \brief A strongly typed enum class representing the vector affordance
*/
  MAKE_ENUM(VectorAffordances,
            noAffordance,
            forces,
            microMagneticForces
      );

  /**
* MatrixAffordances
* \ingroup Affordancestags
* \brief A strongly typed enum class representing the matrix affordance
*/
  MAKE_ENUM(MatrixAffordances,
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

  /**
*
* \ingroup FEParameterTags
* \brief A strongly typed enum class representing the type of the result request
*/
  MAKE_ENUM(ResultType,
            noType,
            magnetization,
            gradientNormOfMagnetization,
            vectorPotential,
            divergenceOfVectorPotential,
            BField,
            HField,
            cauchyStress,
            PK2Stress,
            linearStress,
            director
      );

// clang-format on

/**
 * \brief Struct representing a collection of affordances.
 */
struct AffordanceCollectionImpl
{
  ScalarAffordances scalarAffordances{ScalarAffordances::noAffordance};
  VectorAffordances vectorAffordances{VectorAffordances::noAffordance};
  MatrixAffordances matrixAffordances{MatrixAffordances::noAffordance};
};

/**
 * \brief Concept to check if a given type is one of the predefined affordance enums or the AffordanceCollectionImpl.
 */
template <typename Type>
concept FEAffordance = std::is_same_v<std::remove_cvref_t<Type>, ScalarAffordances> or
                       std::is_same_v<std::remove_cvref_t<Type>, VectorAffordances> or
                       std::is_same_v<std::remove_cvref_t<Type>, MatrixAffordances> or
                       std::is_same_v<std::remove_cvref_t<Type>, AffordanceCollectionImpl>;

inline constexpr VectorAffordances forces = VectorAffordances::forces;

inline constexpr MatrixAffordances stiffness                   = MatrixAffordances::stiffness;
inline constexpr MatrixAffordances stiffnessdiffBucklingVector = MatrixAffordances::stiffnessdiffBucklingVector;
inline constexpr MatrixAffordances mass                        = MatrixAffordances::mass;

inline constexpr ScalarAffordances potentialEnergy = ScalarAffordances::mechanicalPotentialEnergy;

namespace AffordanceCollections {
  inline constexpr AffordanceCollectionImpl elastoStatics = {ScalarAffordances::mechanicalPotentialEnergy,
                                                             VectorAffordances::forces, MatrixAffordances::stiffness};
}

namespace Impl {
  template <typename T>
  struct DeduceRawVectorType
  {
    static_assert(!std::is_same<T, T>::value, "You should end up in the provided specializations");
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
 * \tparam SolutionVectorType_ Type of the solution vector, defaulting to std::reference_wrapper<Eigen::VectorXd>.
 * \tparam ParameterType_ Type of the parameter, defaulting to std::reference_wrapper<double>.
 *
 */
template <typename SolutionVectorType_ = std::reference_wrapper<Eigen::VectorXd>,
          typename ParameterType_      = std::reference_wrapper<double>>
class FERequirements
{
public:
  using SolutionVectorType    = SolutionVectorType_;
  using SolutionVectorTypeRaw = typename Impl::DeduceRawVectorType<std::remove_cvref_t<SolutionVectorType_>>::Type;
  using ParameterType         = ParameterType_;
  using ParameterTypeRaw      = typename ParameterType_::type;

  /**
   * \brief Add an affordance to the requirements.
   *
   * This function adds the specified affordance to the requirements.
   *
   * \tparam Affordance Type of affordance to be added.
   * \param affordance The affordance to be added.
   * \return Reference to the updated FERequirements instance.
   */
  template <FEAffordance Affordance>
  FERequirements& addAffordance(Affordance&& affordance) {
    if constexpr (std::is_same_v<Affordance, ScalarAffordances>)
      affordances.scalarAffordances = affordance;
    else if constexpr (std::is_same_v<Affordance, VectorAffordances>)
      affordances.vectorAffordances = affordance;
    else if constexpr (std::is_same_v<Affordance, MatrixAffordances>)
      affordances.matrixAffordances = affordance;
    else if constexpr (std::is_same_v<Affordance, AffordanceCollectionImpl>)
      affordances = affordance;
    return *this;
  }

  /**
   * \brief Insert a parameter into the requirements.
   *
   * This function inserts the specified parameter into the requirements.
   *
   * \param key The key representing the parameter.
   * \param val Reference to the raw parameter value.
   * \return Reference to the updated FERequirements instance.
   */
  FERequirements& insertParameter(const FEParameter& key, ParameterTypeRaw& val) {
    parameter.insert_or_assign(key, val);
    return *this;
  }

  /**
   * \brief Insert a global solution vector into the requirements.
   *
   * This function inserts the specified global solution vector into the requirements.
   *
   * \param key The key representing the type of the solution vector.
   * \param sol Reference to the raw global solution vector.
   * \return Reference to the updated FERequirements instance.
   */
  FERequirements& insertGlobalSolution(const FESolutions& key, SolutionVectorTypeRaw& sol) {
    sols.insert_or_assign(key, sol);
    return *this;
  }

  /**
   * \brief Get the raw global solution vector for a specific type.
   *
   * This function retrieves the raw global solution vector for the specified type.
   *
   * \param key The key representing the type of the solution vector.
   * \return Const reference to the raw global solution vector.
   *
   * \throws Dune::RangeError if the specified type is not found in the requirements.
   */
  const SolutionVectorTypeRaw& getGlobalSolution(const FESolutions& key) const {
    try {
      if constexpr (std::is_same_v<SolutionVectorType, std::reference_wrapper<Eigen::VectorXd>>)
        return sols.at(key).get();
      else
        return sols.at(key);
    } catch (std::out_of_range& oor) {
      DUNE_THROW(Dune::RangeError, std::string("Out of Range error: ") + std::string(oor.what()) +
                                       " in getGlobalSolution with key" + toString(key));
      abort();
    }
  }

  /**
   * \brief Get the raw parameter value for a specific key.
   *
   * This function retrieves the raw parameter value for the specified key.
   *
   * \param key The key representing the parameter.
   * \return Const reference to the raw parameter value.
   *
   * \throws Dune::RangeError if the specified key is not found in the requirements.
   */
  const ParameterTypeRaw& getParameter(FEParameter&& key) const {
    try {
      return parameter.at(key).get();
    } catch (std::out_of_range& oor) {
      DUNE_THROW(Dune::RangeError, std::string("Out of Range error: ") + std::string(oor.what()) +
                                       " in getParameter with key" + toString(key));
      abort();
    }
  }

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
  bool hasAffordance(Affordance&& affordance) const {
    if constexpr (std::is_same_v<Affordance, ScalarAffordances>)
      return affordances.scalarAffordances == affordance;
    else if constexpr (std::is_same_v<Affordance, VectorAffordances>)
      return affordances.vectorAffordances == affordance;
    else if constexpr (std::is_same_v<Affordance, MatrixAffordances>)
      return affordances.matrixAffordances == affordance;
    else if constexpr (std::is_same_v<Affordance, AffordanceCollectionImpl>)
      return affordances == affordance;
  }

private:
  std::map<FESolutions, SolutionVectorType> sols;
  std::map<FEParameter, ParameterType> parameter;
  AffordanceCollectionImpl affordances;
};

/**
 * \class FErequirements
 * \brief Class representing the requirements for finite element calculations.
 * \deprecated FErequirements is deprecaded and will be removed after v0.5. Use FERequirements instead.
 *
 * This class defines the requirements for finite element calculations, including the types of solution vectors
 * and parameters needed. It provides methods to add affordances, insert parameters, and manage global solution
 * vectors.
 *
 * \tparam SolutionVectorType_ Type of the solution vector, defaulting to std::reference_wrapper<Eigen::VectorXd>.
 * \tparam ParameterType_ Type of the parameter, defaulting to std::reference_wrapper<double>.
 *
 */
template <typename SolutionVectorType_ = std::reference_wrapper<Eigen::VectorXd>,
          typename ParameterType_      = std::reference_wrapper<double>>
class [[deprecated(
    "FErequirements is deprecaded and will be removed after v0.5. Use FERequirements instead.")]] FErequirements
    : public FERequirements<SolutionVectorType_, ParameterType_>
{
  using Base = FERequirements<SolutionVectorType_, ParameterType_>;

public:
  FErequirements() = default;
  FErequirements(Base&& base)
      : Base(std::forward<Base>(base)) {}
  FErequirements(const Base& base)
      : Base(base) {}
  FErequirements& operator=(const Base& base) {
    Base::operator=(base);
    return *this;
  }
  FErequirements& operator=(Base&& base) {
    Base::operator=(std::forward<Base>(base));
    return *this;
  }
};

template <typename Type>
concept ResultTypeConcept = std::is_same_v<Type, ResultType>;

} // namespace Ikarus
