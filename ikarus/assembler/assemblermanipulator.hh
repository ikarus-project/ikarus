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
#include <ikarus/assembler/interface.hh>
#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/dirichletvalues.hh>

namespace Ikarus {

template <template <typename> typename Wrapper, typename Assembler>
requires true
struct ScalarWrapperBase : public Assembler, public ScalarAssembler<Wrapper<Assembler>, typename Assembler::FEContainer,
                                           typename Assembler::DirichletValuesType, typename Assembler::ScalarType>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;
  using ScalarType       = typename Assembler::ScalarType;
  using Interface        = ScalarAssembler<WrappedAssembler, FEC, DV, ScalarType>;
  friend Interface;
  using FunctionType =
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
  std::vector<FunctionType> sfs;

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
    : public Assembler, public VectorAssembler<Wrapper<Assembler>, typename Assembler::FEContainer,
                             typename Assembler::DirichletValuesType, typename Assembler::VectorType>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;

  using VectorType = typename Assembler::VectorType;
  using Interface  = VectorAssembler<Wrapper<Assembler>, typename Assembler::FEContainer,
                                     typename Assembler::DirichletValuesType, typename Assembler::VectorType>;
  friend Interface;
  using FunctionType =
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
  std::vector<FunctionType> vfs;

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
struct MatrixWrapperBase
    : public Assembler, public MatrixAssembler<Wrapper<Assembler>, typename Assembler::FEContainer,
                             typename Assembler::DirichletValuesType, typename Assembler::MatrixType>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;

  using MatrixType = typename Assembler::MatrixType;
  using Interface  = MatrixAssembler<Wrapper<Assembler>, typename Assembler::FEContainer,
                                     typename Assembler::DirichletValuesType, typename Assembler::MatrixType>;
  friend Interface;
  using FunctionType =
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
  std::vector<FunctionType> mfs;

protected:
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
  using ScalarType       = typename Assembler::ScalarType;
  using VectorType       = typename Assembler::VectorType;
  using MatrixType       = typename Assembler::MatrixType;
  using ScalarWrapper    = ScalarWrapperBase<Wrapper, Assembler>;
  using VectorWrapper    = VectorWrapperBase<Wrapper, Assembler>;
  using MatrixWrapper    = MatrixWrapperBase<Wrapper, Assembler>;
  using ScalarInterface  = typename ScalarWrapper::Interface;
  using VectorInterface  = typename VectorWrapper::Interface;
  using MatrixInterface  = typename MatrixWrapper::Interface;

  using MatrixInterface::matrix;
  using ScalarInterface::scalar;
  using VectorInterface::vector;

  using MatrixWrapper::bind;
  using ScalarWrapper::bind;
  using VectorWrapper::bind;

  using ScalarWrapper::getScalarImpl;

  using VectorWrapper::getRawVectorImpl;
  using VectorWrapper::getReducedVectorImpl;
  using VectorWrapper::getVectorImpl;

  using MatrixWrapper::getMatrixImpl;
  using MatrixWrapper::getRawMatrixImpl;
  using MatrixWrapper::getReducedMatrixImpl;
};

template <template <typename> typename Wrapper, typename Assembler>
struct AssemblerWrapperScalarBase : public ScalarWrapperBase<Wrapper, Assembler>
{
  using WrappedAssembler = Wrapper<Assembler>;
  using FEC              = typename Assembler::FEContainer;
  using DV               = typename Assembler::DirichletValuesType;
  using FERequirement    = typename Assembler::FERequirement;

  using ScalarType      = typename Assembler::ScalarType;
  using ScalarInterface = ScalarAssembler<Wrapper<Assembler>, FEC, DV, ScalarType>;
  using ScalarWrapper   = ScalarWrapperBase<Wrapper, Assembler>;
  using ScalarInterface::scalar;
  using ScalarWrapper::bind;
  using ScalarWrapper::getScalarImpl;
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
