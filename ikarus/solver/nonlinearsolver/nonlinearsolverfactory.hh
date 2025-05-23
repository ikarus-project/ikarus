// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverfactory.hh
 * \brief Contains the generic NonlinearSolverFactory class.
 */

#pragma once

#include <utility>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/differentiablefunctionfactory.hh>

namespace Ikarus {

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

    using CorrectionType = typename fTraits::template Range<1>;
    using Domain         = typename fTraits::Domain;
    auto updateF         = [assembler, setting = settings]<typename D, typename C>(D& x, const C& b) {
      if (assembler->dBCOption() == DBCOption::Reduced) {
        setting.updateFunction(x, assembler->createFullVector(b));
      } else
        setting.updateFunction(x, b);
    };
    auto settingsNew = settings.rebindUpdateFunction(std::move(updateF));
    return createNonlinearSolver(std::move(settingsNew), std::move(f));
  }
};
}; // namespace Ikarus