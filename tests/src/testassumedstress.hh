// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testcommon.hh"

#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolverfactory.hh>

template <typename FE>
requires(FE::hasAssumedStress() and FE::strainType == Ikarus::StrainTags::linear)
struct ElementTest<FE>
{
  [[nodiscard]] static auto test() {
    auto easFunctor = [](auto& f, auto& fe, auto& req, auto& affordance) {
      const auto& localView = fe.localView();
      const auto& element   = localView.element();
      constexpr int gridDim = std::remove_cvref_t<decltype(element)>::dimension;

      Dune::TestSuite t("EAS specific test");

      auto subOp = derivative(f);
      std::array<int, (gridDim == 2) ? 1 : 3> parameters;
      if constexpr (gridDim == 2)
        parameters = {5};
      else
        parameters = {18, 24, 30};

      Eigen::VectorXd eigenValues;

      for (auto& numberOfASParameter : parameters) {
        fe.setAssumedStressType(numberOfASParameter);
        auto messageIfFailed = "The numbers of PS parameters are " + std::to_string(numberOfASParameter) + ".";
        t.subTest(checkJacobianOfElement(subOp, req, messageIfFailed));

        t.subTest(checkFEByAutoDiff(f, fe, req, affordance, messageIfFailed));

        decltype(auto) stiffnessMatrix = derivative(subOp)(req);

        auto es = Ikarus::makeIdentitySymEigenSolver<Ikarus::EigenValueSolverType::Spectra>(stiffnessMatrix);
        es.compute();
        eigenValues = es.eigenvalues();

        t.check((eigenValues.array() < 1e-6 * eigenValues.norm()).count() == 3 * gridDim - 3, "Number of eigenvalues")
            << "We always should have " << 3 * gridDim - 3 << " zero eigenvalues, for " << 3 * gridDim - 3
            << " rigid body motions"
            << " but we counted " << (eigenValues.array() < 1e-6 * eigenValues.norm()).count() << "\nEigenValues: \n"
            << eigenValues.transpose() << "The tolerance is " << 1e-6 * eigenValues.norm()
            << "The number of PS parameters is" << numberOfASParameter << std::endl;

        t.check(numberOfASParameter == fe.numberOfInternalVariables())
            << "Number of PS Parameters should be " << numberOfASParameter << "but is "
            << fe.numberOfInternalVariables();

        t.checkThrow([&]() { fe.setAssumedStressType(100); }) << "fe.setAssumedStressType(100) should have failed";

        auto e = calculateScalar(fe, req, Ikarus::ScalarAffordance::mechanicalPotentialEnergy);
        t.check(not isnan(e)) << "Potential should not be NaN";
        t.check(Dune::FloatCmp::ne(e, 0.0)) << "Potential should not be zero";
      }
      return t;
    };
    return easFunctor;
  }
};
