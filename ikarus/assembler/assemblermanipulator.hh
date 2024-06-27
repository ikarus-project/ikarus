// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file assemblermanipulator.hh
 * \brief Defines a decorator for the assemblers that helps to manipulate the assembled quantities
 */

#pragma once

#include <utility>

#include <dune/functions/backends/istlvectorbackend.hh>

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
class AssemblerManipulator : public A, public A::template BaseTemplate<A>
{
protected:
  friend typename A::template BaseTemplate<A>;

public:
  using BaseAssembler            = A;
  using FERequirement            = typename A::FERequirement;
  using AffordanceCollectionType = typename A::AffordanceCollectionType;
  using typename BaseAssembler::MatrixType;
  using typename BaseAssembler::ScalarType;
  using typename BaseAssembler::VectorType;

  using scalarFunction = std::function<void(const BaseAssembler&, const FERequirement&, ScalarAffordance, ScalarType&)>;
  using vectorFunction =
      std::function<void(const BaseAssembler&, const FERequirement&, VectorAffordance, DBCOption, VectorType&)>;
  using matrixFunction =
      std::function<void(const BaseAssembler&, const FERequirement&, MatrixAffordance, DBCOption, MatrixType&)>;

  template <typename... Args>
  requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, AssemblerManipulator>)
  explicit AssemblerManipulator(Args&&... args)
      : baseAssembler(std::forward<Args>(args)...) {}

  /**
   * \brief A helper function to add functions that can be used to manipulate the assembled quantity.
   * \tparam F Type of the function
   * \param f A function that manipulates the assembled quantity.
   */
  template <typename F>
  requires(std::convertible_to<F, scalarFunction> or std::convertible_to<F, vectorFunction> or
           std::convertible_to<F, matrixFunction>)
  void bind(F&& f) {
    if constexpr (std::convertible_to<F, scalarFunction>)
      sfs.emplace_back(std::forward<F>(f));
    else if constexpr (std::convertible_to<F, vectorFunction>)
      vfs.emplace_back(std::forward<F>(f));
    else if constexpr (std::convertible_to<F, matrixFunction>)
      mfs.emplace_back(std::forward<F>(f));
    else
      DUNE_THROW(Dune::IOError, "Function type doesn't meet the requirements.");
  }

private:
  ScalarType& getScalarImpl(const FERequirement& feRequirements, ScalarAffordance affordance) {
    auto& sca = baseAssembler.getScalarImpl(feRequirements, affordance);
    for (const auto sf : sfs)
      sf(baseAssembler, feRequirements, affordance, sca);
    return sca;
  }

  const VectorType& getRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = baseAssembler.getRawVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(baseAssembler, feRequirements, affordance, DBCOption::Raw, vec);
    return vec;
  }

  const VectorType& getVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = baseAssembler.getVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(baseAssembler, feRequirements, affordance, DBCOption::Full, vec);
    return vec;
  }

  const VectorType& getReducedVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = baseAssembler.getReducedVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(baseAssembler, feRequirements, affordance, DBCOption::Reduced, vec);
    return vec;
  }

  const MatrixType& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = baseAssembler.getRawMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(baseAssembler, feRequirements, affordance, DBCOption::Raw, mat);
    return mat;
  }

  const MatrixType& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = baseAssembler.getMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(baseAssembler, feRequirements, affordance, DBCOption::Full, mat);
    return mat;
  }

  const MatrixType& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = baseAssembler.getReducedMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(baseAssembler, feRequirements, affordance, DBCOption::Reduced, mat);
    return mat;
  }

  std::vector<scalarFunction> sfs;
  std::vector<vectorFunction> vfs;
  std::vector<matrixFunction> mfs;

  A baseAssembler;
};
} // namespace Ikarus
