// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file mixin.hh
 * \brief Implementation of the finite element CRTP mixin class
 * \ingroup  finiteelements
 */

#pragma once

#include <dune/common/tuplevector.hh>

#include <ikarus/finiteelements/fetraits.hh>

namespace Ikarus {
/**
 * @brief CRTP mixin class for finite elements with additional skills.
 *
 * This mixin class is designed for finite elements and provides the ability to include additional skills
 * through the CRTP (Curiously Recurring Template Pattern).
 *
 * @tparam PreFE The base finite element class.
 * @tparam Skills A template parameter pack for additional skills to be mixed into the finite element.
 */
template <typename PreFE, template <typename, typename> class... Skills>
struct FEMixin : Skills<PreFE, typename PreFE::template FE<Skills...>>...
{
  /**
   * @brief Constructor for the FEMixin class.
   *
   * @param skillsArgs Variadic arguments for initializing the additional skills.
   */
  explicit FEMixin(typename Skills<PreFE, typename PreFE::template FE<Skills...>>::Pre&&... skillsArgs)
      : Skills<PreFE, typename PreFE::template FE<Skills...>>(
            std::forward<typename Skills<PreFE, typename PreFE::template FE<Skills...>>::Pre>(skillsArgs))... {}

  /**
   * @brief Checks if the mixin class has a specific skill.
   *
   * @tparam Skill The skill to check for.
   * @return True if the skill is present, false otherwise.
   */
  template <template <typename, typename> class Skill>
  static consteval bool hasSkill() {
    return Ikarus::traits::hasType<Skill<PreFE, typename PreFE::template FE<Skills...>>,
                                   std::tuple<Skills<PreFE, typename PreFE::template FE<Skills...>>...>>::value;
  }

  /**
   * @brief Type alias for the supported result types by the mixin.
   */
  using SupportedResultTypes = decltype(std::tuple_cat(
      std::declval<typename Skills<PreFE, typename PreFE::template FE<Skills...>>::SupportedResultTypes>()...));

  using Traits                  = PreFE::Traits;
  using FERequirementType       = typename Traits::FERequirementType;
  using LocalView               = typename Traits::LocalView;
  static constexpr int worldDim = Traits::worlddim;

  /**
   * \brief Calculate the scalar value associated with the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param par The FERequirementType object specifying the requirements for the calculation.
   * \return The calculated scalar value.
   */
  friend auto calculateScalar(const FEMixin& self, const FERequirementType& par) {
    return self.template calculateScalarImpl<double>(par);
  }

  /**
   * \brief Calculate the vector associated with the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param par The FERequirementType object specifying the requirements for the calculation.
   * \param force The vector to store the calculated result.
   */
  friend void calculateVector(const FEMixin& self, const FERequirementType& par,
                              typename Traits::template VectorType<> force) {
    self.template calculateVectorImpl<double>(par, force);
  }

  /**
   * \brief Calculate the matrix associated with the given FERequirementType.
   *
   * \tparam ScalarType The scalar type for the calculation.
   * \param par The FERequirementType object specifying the requirements for the calculation.
   * \param K The matrix to store the calculated result.
   */
  friend void calculateMatrix(const FEMixin& self, const FERequirementType& par,
                              typename Traits::template MatrixType<> K) {
    self.template calculateMatrixImpl<double>(par, K);
  }

  using Skills<PreFE, typename PreFE::template FE<Skills...>>::calculateAtImpl...;

  /**
   * @brief Calculate the element values at a specific location for a given ResultType.
   *
   * @tparam RT The ResultType to calculate.
   * @param req The FERequirementType object specifying the requirements for the calculation.
   * @param local The local coordinates where the calculation is performed.
   * @return The calculated result as specified by the ResultType.
   */
  template <template <typename, int, int> class RT>
  requires requires(FEMixin m, const FERequirementType& req, const Dune::FieldVector<double, Traits::mydim>& local) {
    m.template calculateAtImpl<RT>(req, local, Dune::PriorityTag<10>());
  }
  auto calculateAt(const FERequirementType& req, const Dune::FieldVector<double, Traits::mydim>& local) const {
    return this->template calculateAtImpl<RT>(req, local, Dune::PriorityTag<10>());
  }

  /**
   * @brief Call all bind functions if the skill implement it
   */
  void bind() {
    auto skVisitor = []<typename Skill>(Skill& skill) {
      if constexpr (requires { skill.bind(); })
        skill.bind();
    };
    (skVisitor(static_cast<Skills<PreFE, typename PreFE::template FE<Skills...>>&>(*this)), ...);
  }

  /**
   * @brief Calculate the scalar value in each skill and joins them by `+`.
   *
   * @tparam ScalarType The scalar type for the calculation.
   * @param par The FERequirementType object specifying the requirements for the calculation.
   * @param dx Optional vector used in the calculation.
   * @return The calculated scalar value.
   */
  template <typename ScalarType = double>
  auto calculateScalarImpl(
      const FERequirementType& par,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    return (Skills<PreFE, typename PreFE::template FE<Skills...>>::template calculateScalarImpl<ScalarType>(par, dx) +
            ... + ScalarType{0});
  }

  /**
   * @brief Calculate the vector for each skill
   *
   * @tparam ScalarType The scalar type for the calculation.
   * @param par The FERequirementType object specifying the requirements for the calculation.
   * @param force The vector to store the calculated result.
   * @param dx Optional vector used in the calculation.
   */
  template <typename ScalarType>
  void calculateVectorImpl(
      const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    (Skills<PreFE, typename PreFE::template FE<Skills...>>::template calculateVectorImpl<ScalarType>(par, force, dx),
     ...);
  }

  /**
   * @brief Calculate the matrix for each skill
   *
   * @tparam ScalarType The scalar type for the calculation.
   * @param par The FERequirementType object specifying the requirements for the calculation.
   * @param K The matrix to store the calculated result.
   * @param dx Optional vector used in the calculation.
   */
  template <typename ScalarType>
  void calculateMatrixImpl(
      const FERequirementType& par, typename Traits::template MatrixType<ScalarType> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    (Skills<PreFE, typename PreFE::template FE<Skills...>>::template calculateMatrixImpl<ScalarType>(par, K, dx), ...);
  }

protected:
  /**
   * @brief Get a reference to the underlying finite element object.
   *
   * @return A reference to the underlying finite element object.
   */
  const auto& underlying() const { return static_cast<const typename PreFE::template FE<Skills...>&>(*this); }
  /**
   * @brief Get a reference to the underlying finite element object.
   *
   * @return A reference to the underlying finite element object.
   */
  auto& underlying() { return static_cast<typename PreFE::template FE<Skills...>&>(*this); }
};

/**
 * @brief Struct representing a collection of skills.
 *
 * @tparam ARGS Variadic template parameters representing the skills.
 */
template <typename... ARGS>
struct Skills
{
  using Args = std::tuple<ARGS...>;
  Args args;
};

/**
 * @brief Function to create a Skills instance with the given skills.
 *
 * @tparam Args Variadic template parameters representing the skills.
 * @param args Variadic arguments representing the skills.
 * @return A Skills instance containing the specified skills.
 */
template <typename... Args>
auto skills(Args&&... args) {
  return Skills<Args...>{std::forward_as_tuple(std::decay_t<Args>(args)...)};
}

/**
 * @brief Function to merge two Skills instances.
 *
 * @tparam Args1 Variadic template parameters representing the skills of the first instance.
 * @tparam Args2 Variadic template parameters representing the skills of the second instance.
 * @param sk1 The first Skills instance.
 * @param sk2 The second Skills instance.
 * @return A new Skills instance containing the merged skills.
 */
template <typename... Args1, typename... Args2>
auto merge(Skills<Args1...>&& sk1, Skills<Args2...>&& sk2) {
  return Skills<Args1..., Args2...>{std::tuple_cat(std::move(sk1.args), std::move(sk2.args))};
}

} // namespace Ikarus
