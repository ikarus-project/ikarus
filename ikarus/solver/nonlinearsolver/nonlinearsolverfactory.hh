// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverfactory.hh
 * \brief Contains the generic NonlinearSolverFactory class.
 */

#pragma once

#include <type_traits>
#include <utility>

#include <dune/common/float_cmp.hh>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/differentiablefunctionfactory.hh>
#include <ikarus/utils/functionhelper.hh>

namespace Ikarus {

namespace Impl {
  struct IDBCForceFunction
  {
    template <typename A>
    auto operator()(const A& assembler) {
      return [&]<typename D>(const D&) { return utils::obtainForcesDueToIDBC(assembler); };
    }
  };

  template <typename CT, typename A, typename S>
  auto updateFunctor(const A& assembler, const S& setting) {
    return [&]<typename D, typename C>(D& x, const C& b) {
      if constexpr (not std::is_same_v<C, utils::SyncFERequirements>) {
        // the right-hand side is reduced
        if (assembler->dBCOption() == DBCOption::Reduced and assembler->reducedSize() == b.size()) {
          setting.updateFunction(x, assembler->createFullVector(b));
        } else if (assembler->reducedSize() == x.size() and
                   assembler->size() == b.size()) // the right-hand side is full but x is reduced
        {
          setting.updateFunction(x, assembler->createReducedVector(b));
        } else
          setting.updateFunction(x, b);
      } else { // updates due to inhomogeneous bcs
        if constexpr (requires { x.parameter(); }) {
          auto& dv  = assembler->dirichletValues();
          CT newInc = CT::Zero(dv.size());
          dv.evaluateInhomogeneousBoundaryCondition(newInc, x.parameter());
          for (const auto i : Dune::range(newInc.size()))
            if (Dune::FloatCmp::ne(newInc[i], 0.0))
              x.globalSolution()[i] = newInc[i];
        }
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

  /**
   * \brief A helper function to create another NonlinearSolverFactory object after binding IDBCForceFunction to the
   * configuration of the nonlinear solver. This then helps to handle inhomogeneous Dirichlet BCs when using a certain
   * control routine to perform a nonlinear analysis.
   *
   * \tparam Assembler The type of the assembler used for creating the nonlinear solver.
   *
   * \param assembler The assembler to be used for creating the nonlinear solver.
   *
   * \return The created
   * nonlinear solver factory.
   */
  template <typename Assembler>
  auto withIDBCForceFunction(Assembler&& assembler) const {
    auto idbcForceF  = Impl::IDBCForceFunction{}.template operator()(assembler);
    auto newSettings = settings.rebindIDBCForceFunction(std::move(idbcForceF));
    return NonlinearSolverFactory<std::decay_t<decltype(newSettings)>>{std::move(newSettings)};
  }

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
   *       is used based on whether the option is set to Reduced or not. If inhomogeneous Dirichlet BCs are present,
   *       please call withIDBCForceFunction first.
   */
  template <typename Assembler>
  requires Concepts::FlatAssembler<typename std::remove_cvref_t<Assembler>::element_type>
  auto create(Assembler&& assembler) const {
    auto f        = DifferentiableFunctionFactory::op(assembler);
    using fTraits = typename decltype(f)::Traits;

    using CorrectionType = std::remove_cvref_t<typename fTraits::template Range<1>>;
    using Domain         = typename fTraits::Domain;

    auto updateF     = Impl::updateFunctor<CorrectionType>(assembler, settings);
    auto settingsNew = settings.rebindUpdateFunction(std::move(updateF));
    return createNonlinearSolver(std::move(settingsNew), std::move(f));
  }

private:
  NLSSetting settings;
};
} // namespace Ikarus