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
class AssemblerManipulator : private A, public A::template BaseTemplate<AssemblerManipulator<A>>
{
protected:
  friend typename A::template BaseTemplate<AssemblerManipulator>;

  using Base = A::template BaseTemplate<AssemblerManipulator<A>>;

public:
  using Assembler            = A;
  using FERequirement            = typename A::FERequirement;
  using AffordanceCollectionType = typename A::AffordanceCollectionType;
  using typename Assembler::MatrixType;
  using typename Assembler::ScalarType;
  using typename Assembler::VectorType;

  using Base::matrix;
  using Base::vector;
  using Base::scalar;

  using scalarFunction = std::function<void(const AssemblerManipulator&, const FERequirement&, ScalarAffordance, ScalarType&)>;
  using vectorFunction =
      std::function<void(const AssemblerManipulator&, const FERequirement&, VectorAffordance, DBCOption, VectorType&)>;
  using matrixFunction =
      std::function<void(const AssemblerManipulator&, const FERequirement&, MatrixAffordance, DBCOption, MatrixType&)>;

  template <typename... Args>
  requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, AssemblerManipulator>)
  explicit AssemblerManipulator(Args&&... args)
      : Assembler(std::forward<Args>(args)...),Base(std::forward<Args>(args)...) {}

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
    auto& sca = Assembler::getScalarImpl(feRequirements, affordance);
    for (const auto sf : sfs)
      sf(*this, feRequirements, affordance, sca);
    return sca;
  }

  const VectorType& getRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = Assembler::getRawVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(*this, feRequirements, affordance, DBCOption::Raw, vec);
    return vec;
  }

  const VectorType& getVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = Assembler::getVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(*this, feRequirements, affordance, DBCOption::Full, vec);
    return vec;
  }

  const VectorType& getReducedVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance) {
    auto& vec = Assembler::getReducedVectorImpl(feRequirements, affordance);
    for (const auto vf : vfs)
      vf(*this, feRequirements, affordance, DBCOption::Reduced, vec);
    return vec;
  }

  const MatrixType& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = Assembler::getRawMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(*this, feRequirements, affordance, DBCOption::Raw, mat);
    return mat;
  }

  const MatrixType& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = Assembler::getMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(*this, feRequirements, affordance, DBCOption::Full, mat);
    return mat;
  }

  const MatrixType& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance) {
    MatrixType& mat = Assembler::getReducedMatrixImpl(feRequirements, affordance);
    for (const auto mf : mfs)
      mf(*this, feRequirements, affordance, DBCOption::Reduced, mat);
    return mat;
  }

  std::vector<scalarFunction> sfs;
  std::vector<vectorFunction> vfs;
  std::vector<matrixFunction> mfs;

  // A baseAssembler;
};
} // namespace Ikarus
