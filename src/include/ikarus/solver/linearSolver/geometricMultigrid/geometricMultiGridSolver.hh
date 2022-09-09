//
// Created by ac129893 on 30.06.2022.
//

#pragma once

#include "gridTransfer.h"
#include "src/include/ikarus/finiteElements/feRequirements.hh"
#include "src/include/ikarus/solver/linearSolver/linearSolver.hh"

#include <iostream>
#include <utility>

#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <Eigen/Core>

#include "ikarus/assembler/simpleAssemblers.hh"

namespace Ikarus {
  template <typename Grid, typename PreBasisFactory, typename FEContainer>
  class GeometricMultiGridSolver {
    using RequirementType = typename FEContainer::value_type::FERequirementType;

    using BasisType
        = decltype(makeBasis(std::declval<typename Grid::LevelGridView>(), std::declval<PreBasisFactory>()));

  public:
    GeometricMultiGridSolver(const Grid* grid, const PreBasisFactory& preBasisFactory,
                             const FEContainer& feVectorCoarse)
        : transfer{grid}, directSolver{SolverTypeTag::sd_SimplicialLDLT} {
      fes.resize(grid->maxLevel() + 1);
      finestLevel = grid->maxLevel();
      fes[0]      = feVectorCoarse;
      for (int level = 1; level < grid->maxLevel() + 1; ++level) {
        auto coarseGridView = grid->levelGridView(level - 1);
        auto fineGridView   = grid->levelGridView(level);

        auto coarseBasis = makeBasis(coarseGridView, preBasisFactory);
        auto fineBasis   = makeBasis(fineGridView, preBasisFactory);

        auto coarseElement = elements(coarseGridView).begin();
        fes[level].reserve(fineGridView.size(0));
        for (auto& coarseFe : fes[level - 1]) {
          {
            for (auto& childsElement : descendantElements(*coarseElement, coarseElement->level() + 1)) {
              fes[level].emplace_back(fineBasis, childsElement, coarseFe.settings());
            }
            ++coarseElement;
          }
        }
      }

      for (int level = 0; level < grid->maxLevel() + 1; ++level) {
        auto gridView = grid->levelGridView(level);

        auto basis          = std::make_shared<BasisType>(makeBasis(gridView, preBasisFactory));
        auto dirichletFlags = std::make_shared<std::vector<bool>>(basis->size(), false);

        Dune::Functions::forEachBoundaryDOF(*basis, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
          if (std::abs(intersection.geometry().center()[1]) < 1e-8) {
            (*dirichletFlags)[localView.index(localIndex)[0]] = true;
          }
        });
        assemblers.emplace_back(basis, fes[level], dirichletFlags);
      }

      transfer.createOperators(preBasisFactory);
      iterativeSolver.setMaxIterations(1);
    }

    mutable Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                                     Eigen::IncompleteCholesky<double>>
        iterativeSolver;
    Ikarus::GridTransfer<Grid> transfer;
    mutable ILinearSolver<double> directSolver;
    mutable RequirementType requirementType;

    mutable std::vector<Ikarus::SparseFlatAssembler<BasisType, FEContainer>> assemblers;
    std::vector<FEContainer> fes;
    std::vector<BasisType> bases;
    int finestLevel{0};

  private:
    template <typename ScalarType>
    void smoothing(const Eigen::VectorX<ScalarType>& dFineFull, const Eigen::VectorX<ScalarType>& RfineRed,
                   Eigen::VectorX<ScalarType>& dFineRed) const {
      assemblers[finestLevel].createReducedVector(dFineFull, dFineRed);
      dFineRed = iterativeSolver.solveWithGuess(RfineRed, dFineRed);
    }

  public:
    void transformToFineFull(const Eigen::VectorX<double>& dFineRed, Eigen::VectorX<double>& dFineFull) {
      assemblers[finestLevel].createFullVector(dFineRed, dFineFull);
    }

    template <typename ScalarType>
    void solve(Eigen::VectorX<ScalarType>& dFineRed, const Eigen::VectorX<ScalarType>& RfineRed_) const {
      //    Eigen::MatrixXd KcoarseRed,KfineRed;
      Eigen::VectorXd RcoarseRed, dCoarseRed, dCoarseFull, RfineFull, RcoarseFull;
      dCoarseFull.setZero(assemblers[finestLevel - 1].size());
      Eigen::VectorXd dFineFull;
      dFineFull.setZero(assemblers[finestLevel].size());

      const double lambdaLoad = -1;

      requirementType = Ikarus::FErequirementsBuilder()
                            .insertGlobalSolution(Ikarus::FESolutions::displacement, dFineFull)
                            .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                            .addAffordance(Ikarus::VectorAffordances::forces)
                            .build();
      const Eigen::VectorXd& RfineRed = -assemblers[finestLevel].getReducedVector(requirementType);
      requirementType          = Ikarus::FErequirementsBuilder()
                            .insertGlobalSolution(Ikarus::FESolutions::displacement, dCoarseFull)
                            .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                            .addAffordance(Ikarus::MatrixAffordances::stiffness)
                            .build();
      auto& KcoarseRed = assemblers[finestLevel - 1].getReducedMatrix(requirementType);

      directSolver.compute(KcoarseRed);
      assemblers[finestLevel].createFullVector(RfineRed, RfineFull);
      transfer.restrictTo(finestLevel - 1, RfineFull, RcoarseFull);
      assemblers[finestLevel - 1].createReducedVector(RcoarseFull, RcoarseRed);

      directSolver.solve(dCoarseRed, RcoarseRed);

      assemblers[finestLevel - 1].createFullVector(dCoarseRed, dCoarseFull);

      requirementType = Ikarus::FErequirementsBuilder()
                            .insertGlobalSolution(Ikarus::FESolutions::displacement, dFineFull)
                            .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                            .addAffordance(Ikarus::MatrixAffordances::stiffness)
                            .build();
      const auto KfineRed = assemblers[finestLevel].getReducedMatrix(requirementType);
      std::cout<<"BeginKfineRed.template block(2,2,0,0).eval()"<<std::endl;
      std::cout<<KfineRed. block(0,0,2,2)<<std::endl;
      transfer.prolongateFrom(finestLevel - 1, dCoarseFull, dFineFull);
//      std::cout<<KfineRed<<std::endl;
//      std::cout<<"RfineRed.transpose()"<<RfineRed.transpose()<<std::endl;
      iterativeSolver.compute(KfineRed);
      assemblers[finestLevel].createReducedVector(dFineFull, dFineRed);

      Eigen::VectorXd eFineFull, residualMGFineRed, residualMGFineFull, residualMGCoarseFull, eCoarseRed,
          residualMGCoarseRed, eCoarseFull;
      eFineFull.resizeLike(dFineFull);
      residualMGFineRed.resizeLike(dFineRed);
      eFineFull.setOnes();
      residualMGFineRed.setOnes();
//
      int iter          = 0;
      int maxIterations = 1000;
      spdlog::info("iter ResidualNorm: CorrectionNorm");
      while (residualMGFineRed.norm() > 1e-11) {
        smoothing(dFineFull, RfineRed, dFineRed);  // Pre-Smoothing

        residualMGFineRed = RfineRed - KfineRed * dFineRed;
        assemblers[finestLevel].createFullVector(residualMGFineRed, residualMGFineFull);

        transfer.restrictTo(finestLevel - 1, residualMGFineFull, residualMGCoarseFull);

        assemblers[finestLevel - 1].createReducedVector(residualMGCoarseFull, residualMGCoarseRed);
        directSolver.solve(eCoarseRed, residualMGCoarseRed);
        assemblers[finestLevel - 1].createFullVector(eCoarseRed, eCoarseFull);

        transfer.prolongateFrom(finestLevel - 1, eCoarseFull, eFineFull);

        assemblers[finestLevel].createFullVector(dFineRed, dFineFull);
        dFineFull += eFineFull;

        smoothing(dFineFull, RfineRed, dFineRed);  // Post smoothing
        assemblers[finestLevel].createFullVector(dFineRed, dFineFull);
        spdlog::info("{:>6d} {:>9.2e} {:>9.2e}", iter, residualMGFineRed.norm(), eFineFull.norm());

        ++iter;
      }
      dFineFull.setZero();
////
//      requirementType = Ikarus::FErequirementsBuilder()
//          .insertGlobalSolution(Ikarus::FESolutions::displacement, dFineFull)
//          .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
//          .addAffordance(Ikarus::MatrixAffordances::stiffness)
//          .build();
//      std::cout<<"BAKfineRed.template block(2,2,0,0).eval()"<<std::endl;
//      std::cout<<KfineRed. block(0,0,2,2)<<std::endl;
      const auto& KfineRed2 = assemblers[finestLevel].getReducedMatrix(requirementType);
//
      std::cout<<"KfineRed.template block(2,2,0,0).eval()"<<std::endl;
      std::cout<<KfineRed. block(0,0,2,2)<<std::endl;
      std::cout<<"KfineRed2.template block(2,2,0,0).eval()"<<std::endl;
      std::cout<<KfineRed2. block(0,0,2,2)<<std::endl;
      requirementType = Ikarus::FErequirementsBuilder()
          .insertGlobalSolution(Ikarus::FESolutions::displacement, dFineFull)
          .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
          .addAffordance(Ikarus::VectorAffordances::forces)
          .build();
      const Eigen::VectorXd& RfineRed2 = -assemblers[finestLevel].getReducedVector(requirementType);
            Eigen::SimplicialLDLT<std::remove_cvref_t<decltype(KfineRed2)>, Eigen::Lower | Eigen::Upper> solver;
            solver.compute(KfineRed);
            dFineRed = solver.solve(RfineRed);
            assemblers[finestLevel].createFullVector(dFineRed, dFineFull);
    }
  };
}  // namespace Ikarus
