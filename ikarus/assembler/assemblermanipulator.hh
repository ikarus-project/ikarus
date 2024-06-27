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
class AssemblerManipulator : public A::Base
{
protected:
  friend typename A::Base;

public:
  using Assembler                = A;
  using AssemblerBase            = typename A::Base;
  using FERequirement            = typename A::FERequirement;
  using AffordanceCollectionType = typename A::AffordanceCollectionType;
  using typename AssemblerBase::MatrixType;
  using typename AssemblerBase::ScalarType;
  using typename AssemblerBase::VectorType;

  using scalarFunction = std::function<void(const Assembler&, const FERequirement&, ScalarAffordance, ScalarType&)>;
  using vectorFunction =
      std::function<void(const Assembler&, const FERequirement&, VectorAffordance, DBCOption, VectorType&)>;
  using matrixFunction =
      std::function<void(const Assembler&, const FERequirement&, MatrixAffordance, DBCOption, MatrixType&)>;

  template <typename... Args>
  requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, AssemblerManipulator>)
  explicit AssemblerManipulator(Args&&... args)
      : assembler(std::forward<Args>(args)...) {}

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
    auto& sca = assembler.getScalarImpl(feRequirements, affordance);
    for (const auto sf : sfs)
      sf(assembler, feRequirements, affordance, sca);
    return sca;
  }

  const VectorType& getRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = assembler.getRawVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(assembler, feRequirements, affordance, DBCOption::Raw, vec);
    return vec;
  }

  const VectorType& getVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = assembler.getVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(assembler, feRequirements, affordance, DBCOption::Full, vec);
    return vec;
  }

  const VectorType& getReducedVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = assembler.getReducedVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(assembler, feRequirements, affordance, DBCOption::Reduced, vec);
    return vec;
  }

  const MatrixType& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = assembler.getRawMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(assembler, feRequirements, affordance, DBCOption::Raw, mat);
    return mat;
  }

  const MatrixType& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = assembler.getMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(assembler, feRequirements, affordance, DBCOption::Full, mat);
    return mat;
  }

  const MatrixType& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = assembler.getReducedMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(assembler, feRequirements, affordance, DBCOption::Reduced, mat);
    return mat;
  }

  std::vector<scalarFunction> sfs;
  std::vector<vectorFunction> vfs;
  std::vector<matrixFunction> mfs;

  A assembler;
};
} // namespace Ikarus
