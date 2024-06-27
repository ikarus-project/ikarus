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
#include <ikarus/assembler/simpleassemblers.hh> //TODO ALEX TARUN chagne this to interface.hh
#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/dirichletvalues.hh>

namespace Ikarus {

template <template <typename> typename Wrapper, typename Assembler>
requires true
struct ScalarWrapperBase
    : Assembler,
      public ScalarAssemblerBase<Wrapper<Assembler>, typename Assembler::FEContainer,
                                 typename Assembler::DirichletValuesType, typename Assembler::ScalarType>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;

  using ScalarType      = Assembler::ScalarType;
  using ScalarInterface = ScalarAssemblerBase<Wrapper<Assembler>, FEC, DV, ScalarType>;
  friend ScalarInterface;
  using scalarFunction =
      std::function<void(const Wrapper<Assembler>&, const FERequirement&, ScalarAffordance, ScalarType&)>;
  /**
   * \brief A helper function to add functions that can be used to manipulate the assembled quantity.
   * \tparam F Type of the function
   * \param f A function that manipulates the assembled quantity.
   */
  template <typename F>
  requires Concepts::IsFunctorWithArgs<F, const Wrapper<Assembler>&, const FERequirement&, ScalarAffordance,
                                       ScalarType&>
  void bind(F&& f) {
    sfs.emplace_back(std::forward<F>(f));
  }
  std::vector<scalarFunction> sfs;

protected:
  ScalarType& getScalarImpl(const FERequirement& feRequirements, ScalarAffordance affordance) {
    auto& sca = Assembler::getScalarImpl(feRequirements, affordance);
    for (const auto sf : sfs)
      sf(*this, feRequirements, affordance, sca);
    return sca;
  }
};

template <template <typename> typename Wrapper, typename Assembler>
requires true
struct VectorWrapperBase
    : Assembler,
      public VectorAssemblerBase<Wrapper<Assembler>, typename Assembler::FEContainer,
                                 typename Assembler::DirichletValuesType, typename Assembler::VectorType>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;

  using VectorType      = Assembler::VectorType;
  using VectorInterface = VectorAssemblerBase<Wrapper<Assembler>, typename Assembler::FEContainer,
                                              typename Assembler::DirichletValuesType, typename Assembler::VectorType>;
  friend VectorInterface;
  using vectorFunction =
      std::function<void(const Wrapper<Assembler>&, const FERequirement&, VectorAffordance, DBCOption, VectorType&)>;
  /**
   * \brief A helper function to add functions that can be used to manipulate the assembled quantity.
   * \tparam F Type of the function
   * \param f A function that manipulates the assembled quantity.
   */
  template <typename F>
  requires Concepts::IsFunctorWithArgs<F, const Wrapper<Assembler>&, const FERequirement&, VectorAffordance, DBCOption,
                                       VectorType&>
  void bind(F&& f) {
    vfs.emplace_back(std::forward<F>(f));
  }
  std::vector<vectorFunction> vfs;

protected:
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
};

template <template <typename> typename Wrapper, typename Assembler>
requires true
struct MatrixWrapperBase : public Assembler,
                           MatrixAssemblerBase<Wrapper<Assembler>, typename Assembler::FEContainer,
                                               typename Assembler::DirichletValuesType, typename Assembler::MatrixType>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;

  using MatrixType      = Assembler::MatrixType;
  using MatrixInterface = MatrixAssemblerBase<Wrapper<Assembler>, typename Assembler::FEContainer,
                                              typename Assembler::DirichletValuesType, typename Assembler::MatrixType>;
  friend MatrixInterface;
  using matrixFunction =
      std::function<void(const Wrapper<Assembler>&, const FERequirement&, MatrixAffordance, DBCOption, MatrixType&)>;
  /**
   * \brief A helper function to add functions that can be used to manipulate the assembled quantity.
   * \tparam F Type of the function
   * \param f A function that manipulates the assembled quantity.
   */
  template <typename F>
  requires Concepts::IsFunctorWithArgs<F, const Wrapper<Assembler>&, const FERequirement&, MatrixAffordance, DBCOption,
                                       MatrixType&>
  void bind(F&& f) {
    mfs.emplace_back(std::forward<F>(f));
  }
  std::vector<matrixFunction> mfs;

private:
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
};

template <template <typename> typename Wrapper, typename Assembler>
struct AssemblerWrapperAllBase : public ScalarWrapperBase<Wrapper, Assembler>,
                                 public VectorWrapperBase<Wrapper, Assembler>,
                                 public MatrixWrapperBase<Wrapper, Assembler>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;
  using MatrixType       = Assembler::MatrixType;
  using VectorType       = Assembler::VectorType;
  using ScalarType       = Assembler::ScalarType;
  using MatrixInterface  = MatrixAssemblerBase<Wrapper<Assembler>, FEC, DV, MatrixType>;
  using VectorInterface  = VectorAssemblerBase<Wrapper<Assembler>, FEC, DV, VectorType>;
  using ScalarInterface  = ScalarAssemblerBase<Wrapper<Assembler>, FEC, DV, ScalarType>;

  using MatrixInterface::matrix;
  using ScalarInterface::scalar;
  using VectorInterface::vector;

  using MatrixWrapperBase<Wrapper, Assembler>::bind;
  using VectorWrapperBase<Wrapper, Assembler>::bind;
  using ScalarWrapperBase<Wrapper, Assembler>::bind;
};

template <template <typename> typename Wrapper, typename Assembler>
struct AssemblerWrapperScalarBase : private Assembler, public ScalarWrapperBase<Wrapper, Assembler>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;

  using ScalarType      = Assembler::ScalarType;
  using ScalarInterface = ScalarAssemblerBase<Wrapper<Assembler>, FEC, DV, ScalarType>;
  using ScalarInterface::scalar;
};
struct InheritanceDecider
{
  template <template <typename> typename Wrapper, typename Assembler>
  using Type = std::conditional_t<Concepts::FlatAssembler<Assembler>, AssemblerWrapperAllBase<Wrapper, Assembler>,
                                  AssemblerWrapperAllBase<Wrapper, Assembler>>;
};
/**
 * \class AssemblerManipulator
 * \brief The AssemblerManipulator defines a decorator for the assemblers that
 * helps to manipulate the assembled quantities.
 * \ingroup assembler
 *
 * \tparam A Type of the assembler.
 */
template <Concepts::FlatAssembler A>
class AssemblerManipulator : public InheritanceDecider::Type<AssemblerManipulator, A>
{
public:
  using Base = InheritanceDecider::Type<AssemblerManipulator, A>;

  template <typename... Args>
  requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>, AssemblerManipulator>)
  explicit AssemblerManipulator(Args&&... args)
      : Base(std::forward<Args>(args)...) {}
};
} // namespace Ikarus
