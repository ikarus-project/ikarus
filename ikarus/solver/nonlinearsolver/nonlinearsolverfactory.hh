// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverfactory.hh
 * \brief Contains the generic NonlinearSolverFactory class.
 */

#pragma once

#include <type_traits>
#include <utility>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/differentiablefunctionfactory.hh>

namespace Ikarus {

namespace Impl {
  template <typename CT, typename A>
  auto idbcForceFunctor(const A& assembler) {
    return [&]() {
      auto& loadFactor = assembler->requirement().parameter();
      auto& x          = assembler->requirement().globalSolution();
      const auto& K    = assembler->matrix(DBCOption::Raw);
      auto& dv         = assembler->dirichletValues();
      CT newInc        = CT::Zero(dv.size());
      dv.evaluateInhomogeneousBoundaryCondition(newInc, loadFactor);
      dv.setZeroAtConstrainedDofs(x);
      auto F_dirichlet = (K * newInc).eval();
      if (assembler->dBCOption() == DBCOption::Full)
        assembler->dirichletValues().setZeroAtConstrainedDofs(F_dirichlet);
      else
        assembler->createReducedVector(F_dirichlet);

      return F_dirichlet;
    };
  }

  template <typename CT, typename A, typename S>
  auto updateFunctor(const A& assembler, const S& setting) {
    return [&]<typename D, typename C>(D& x, const C& b) {
      if constexpr (not std::is_same_v<C, utils::ZeroIncrementTag>) {
        // the right-hand side is reduced
        if (assembler->dBCOption() == DBCOption::Reduced and assembler->reducedSize() == b.size()) {
          setting.updateFunction(x, assembler->createFullVector(b));
        } else if (assembler->reducedSize() == x.size() and
                   assembler->size() == b.size()) // the right-hand side is full but x is reduced
        {
          setting.updateFunction(x, assembler->createReducedVector(b));
        } else
          setting.updateFunction(x, b);
      }
      // updates due to inhomogeneous bcs
      if constexpr (requires { x.parameter(); }) {
        auto& dv  = assembler->dirichletValues();
        CT newInc = CT::Zero(dv.size());
        dv.evaluateInhomogeneousBoundaryCondition(newInc, x.parameter());
        setting.updateFunction(x, newInc);
      }
    };
  }
} // namespace Impl

/**
 * \brief A factory class for creating nonlinear solvers.
 *
 * This class is responsible for creating nonlinear solvers using the provided settings
 * and an assembler that satisfies the FlatAssembler concept.
 *
 * \tparam NLSSetting The type of the settings used for creating nonlinear solvers.
 */
template <typename NLSSetting>
struct NonlinearSolverFactory
{
  /**
   * \brief Constructs a NonlinearSolverFactory with the given settings.
   *
   * \param s The settings to be used by the factory.
   */
  NonlinearSolverFactory(NLSSetting s)
      : settings(s) {}

  NLSSetting settings;

  /**
   * \brief Creates a nonlinear solver using the provided assembler.
   *
   * The assembler must satisfy the FlatAssembler concept.
   *
   * \tparam Assembler The type of the assembler used for creating the nonlinear solver.
   * \param assembler The assembler to be used for creating the nonlinear solver.
   * \return The created nonlinear solver.
   *
   * \note The assembler's dBCOption is checked, and the appropriate update function
   *       is used based on whether the option is set to Reduced or not.
   */
  template <typename Assembler>
  requires Concepts::FlatAssembler<typename std::remove_cvref_t<Assembler>::element_type>
  auto create(Assembler&& assembler) const {
    auto f        = DifferentiableFunctionFactory::op(assembler);
    using fTraits = typename decltype(f)::Traits;

    using CorrectionType = std::remove_cvref_t<typename fTraits::template Range<1>>;
    using Domain         = typename fTraits::Domain;

    auto updateF        = Impl::updateFunctor<CorrectionType>(assembler, settings);
    auto idbcForceF     = Impl::idbcForceFunctor<CorrectionType>(assembler);
    auto settingsNew    = settings.rebindUpdateFunction(std::move(updateF));
    auto settingsNewNew = settingsNew.rebindIDBCForceFunction(std::move(idbcForceF));
    return createNonlinearSolver(std::move(settingsNewNew), std::move(f));
  }
};
}; // namespace Ikarus