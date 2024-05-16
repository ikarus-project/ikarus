// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <utility>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/nonlinopfactory.hh>

namespace Ikarus
{



template <typename NLSSetting>
struct NonlinearSolverFactory
{
  NonlinearSolverFactory(NLSSetting s)
      : settings(s) {}

  NLSSetting settings;

public:
  template <typename Assembler>
  auto create(Assembler&& assembler) const {
    auto nonLinOp             = NonLinearOperatorFactory::op(assembler);
    std::function updateF     = [assembler, setting = settings](decltype(nonLinOp.firstParameter())& a,
                                                            const decltype(nonLinOp.derivative())& b) {
      if (assembler->enforcingDBCOption() == EnforcingDBCOption::Reduced) {
        setting.updateFunction(a, assembler->createFullVector(b));
      } else
        setting.updateFunction(a, b);
    };
    auto settingsNew = settings.rebindUpdateFunction(std::move(updateF));
    return createNonlinearSolver(std::move(settingsNew), std::move(nonLinOp));
  }
  };
};