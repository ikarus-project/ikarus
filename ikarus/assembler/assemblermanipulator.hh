// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file assemblermanipulator.hh
 * \brief Defines a decorator for the assemblers that helps to manipulate the assembled quantities
 */

#pragma once

#include "simpleassemblers.hh"

#include <utility>

#include <dune/functions/backends/istlvectorbackend.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/dirichletvalues.hh>

namespace Ikarus {
/**
 * \class AssemblerManipulator
 * \brief The AssemblerManipulator defines a decorator for the assemblers that
 * helps to manipulate the assembled quantities.
 * \ingroup assembler
 *
 * \tparam A Type of the assembler.
 */
template <Concepts::FlatAssembler A>
class AssemblerManipulator
{
public:
  using Assembler                = A;
  using FERequirement            = typename A::FERequirement;
  using AffordanceCollectionType = typename A::AffordanceCollectionType;
  using MatrixType =
      std::conditional_t<std::is_same_v<Assembler, SparseFlatAssembler<typename Assembler::FEContainer,
                                                                       typename Assembler::DirichletValuesType>>,
                         Eigen::SparseMatrix<double>, Eigen::MatrixXd>;

  using scalarFunction = std::function<void(const FERequirement&, const ScalarAffordance&, const double&)>;
  using vectorFunction = std::function<void(const FERequirement&, const VectorAffordance&, const Eigen::VectorXd&)>;
  using matrixFunction = std::function<void(const FERequirement&, const MatrixAffordance&, const MatrixType&)>;

  template <typename... Args>
  requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, AssemblerManipulator>)
  explicit AssemblerManipulator(Args&&... args)
      : baseAssembler(std::forward<Args>(args)...) {}

  /**
   * \brief Calculates the scalar quantity requested by feRequirements and affordance.
   *
   * \param feRequirements Reference to the finite element requirements.
   * \param affordance The scalar affordance
   * \return Const reference to the calculated scalar quantity.
   */
  const double& scalar(const FERequirement& feRequirements, ScalarAffordance affordance) {
    double s = baseAssembler.scalar(feRequirements, affordance);
    for (const auto sf : sfs)
      sf(feRequirements, affordance, s);
    return s;
  }

  /**
   * \brief Calculates the scalar quantity requested by the bound feRequirements and returns a reference.
   *
   * \return Const reference to the calculated scalar quantity.
   */
  const double& scalar() {
    return scalar(baseAssembler.requirement(), baseAssembler.affordanceCollection().scalarAffordance());
  }

  /**
    * \brief Calculates the vectorial quantity requested by the  feRequirements and the affordance.
    Depending on the requested DBCOption, the raw, reduced or full vector is returned.
    *
    * \param feRequirements Reference to the finite element requirements.
    * \param affordance The vector affordance
    * \param dbcOption The DBCOption
    * \return Const reference to the calculated vectorial quantity.
    */
  const Eigen::VectorXd& vector(const FERequirement& feRequirements, VectorAffordance affordance,
                                DBCOption dbcOption = DBCOption::Full) {
    const auto& v = baseAssembler.vector(feRequirements, affordance, dbcOption);
    for (const auto vf : vfs)
      vf(feRequirements, affordance, v);
    return v;
  }

  /**
  * \brief Calculates the vectorial quantity requested by the bound feRequirements and the affordance.
  Depending on the requested DBCOption, the raw, reduced or full vector is returned.
  * \param dbcOption The DBCOption
  * \return Const reference to the calculated vectorial quantity.
  */
  const Eigen::VectorXd& vector(DBCOption dbcOption) {
    return vector(baseAssembler.requirement(), baseAssembler.affordanceCollection().vectorAffordance(), dbcOption);
  }

  /**
   * \brief Calculates the vectorial quantity requested by the bound feRequirements, the affordance and the dBCOption.
   * Depending on the DBCOption, the raw, reduced or full vector is returned.
   * \return Const reference to the calculated vectorial quantity.
   */
  const Eigen::VectorXd& vector() { return vector(baseAssembler.dBCOption()); }

  /**
    * \brief  Calculates the matrix quantity requested by feRequirements and the affordance.
    *
    * \param feRequirements Reference to the finite element requirements.
    * \param affordance The matrix affordance
    * \param dbcOption The DBCOption

    * \return Reference to the raw dense matrix quantity.
    */
  const MatrixType& matrix(const FERequirement& feRequirements, MatrixAffordance affordance,
                           DBCOption dbcOption = DBCOption::Full) {
    const auto& m = baseAssembler.matrix(feRequirements, affordance, dbcOption);
    for (const auto mf : mfs)
      mf(feRequirements, affordance, m);
    return m;
  }

  /**
   * \brief  Calculates the matrix quantity requested by the bound  feRequirements and the affordance.
   *
   * \param dbcOption The DBCOption
   * \return Reference to the raw dense matrix quantity.
   */
  const Eigen::MatrixXd& matrix(DBCOption dbcOption) {
    return matrix(baseAssembler.requirement(), baseAssembler.affordanceCollection().matrixAffordance(), dbcOption);
  }

  /**
   * \brief  Calculates the matrix quantity requested by the bound  feRequirements, the affordance and the dBCOption.
   * \return Reference to the dense matrix quantity.
   */
  const Eigen::MatrixXd& matrix() { return matrix(baseAssembler.dBCOption()); }

  /**
   * \brief A helper function to add functions that can be used to manipulate the assembled quantity.
   * \tparam F Type of the function
   * \param f A function that manipulates the assembled quantity.
   */
  template <typename F>
  requires(std::is_same_v<F, scalarFunction> or std::is_same_v<F, vectorFunction> or std::is_same_v<F, matrixFunction>)
  void bind(F&& f) {
    if constexpr (std::is_same_v<F, scalarFunction>)
      sfs.emplace_back(std::forward<F>(f));
    else if constexpr (std::is_same_v<F, vectorFunction>)
      vfs.emplace_back(std::forward<F>(f));
    else
      mfs.emplace_back(std::forward<F>(f));
  }

private:
  std::vector<scalarFunction> sfs;
  std::vector<vectorFunction> vfs;
  std::vector<matrixFunction> mfs;
  A baseAssembler;
};
} // namespace Ikarus
