//
// Created by ac129893 on 30.06.2022.
//

#pragma once

#include "gridTransfer.h"
#include "src/include/ikarus/finiteElements/feRequirements.hh"
#include "src/include/ikarus/solver/linearSolver/linearSolver.hh"

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
      fes.resize(grid->maxLevel()+1);
      finestLevel = grid->maxLevel();
      fes[0] = feVectorCoarse;
      for (int level = 1; level < grid->maxLevel()+1; ++level) {
        auto coarseGridView = grid->levelGridView(level-1);
        auto fineGridView   = grid->levelGridView(level);

        auto coarseBasis = makeBasis(coarseGridView, preBasisFactory);
        auto fineBasis   = makeBasis(fineGridView, preBasisFactory);

        auto coarseElement = elements(coarseGridView).begin();
        fes[level].reserve(fineGridView.size(0));
        for (auto& coarseFe : fes[level-1]) {
          {
            for (auto& childsElement : descendantElements(*coarseElement, coarseElement->level()+1)) {
              fes[level].emplace_back(fineBasis, childsElement, coarseFe.settings());
            }
            ++coarseElement;
          }
        }
      }

      for (int level = 0; level < grid->maxLevel()+1; ++level) {
        auto gridView = grid->levelGridView(level);

        auto basis          = std::make_shared<BasisType>(makeBasis(gridView, preBasisFactory));
        auto dirichletFlags = std::make_shared<std::vector<bool>>((*basis).size(), false);

        Dune::Functions::forEachBoundaryDOF(*basis, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
          if (std::abs(intersection.geometry().center()[1]) < 1e-8) {
            (*dirichletFlags)[localView.index(localIndex)[0]] = true;
          }
        });
        assemblers.emplace_back(basis, fes[level], dirichletFlags);
      }

      transfer.createOperators(preBasisFactory);
      iterativeSolver.setMaxIterations(10000);
    }

    mutable Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> iterativeSolver;
    Ikarus::GridTransfer<Grid> transfer;
    mutable ILinearSolver<double> directSolver;
    mutable RequirementType requirementType;


    mutable std::vector<Ikarus::SparseFlatAssembler<BasisType, FEContainer>> assemblers;
    std::vector<FEContainer> fes;
    std::vector<BasisType> bases;
    int finestLevel{0};

  private:
  public:
    template <typename ScalarType>
    void smoothing(const Eigen::VectorX<ScalarType>& dFineFull, const Eigen::VectorX<ScalarType>& RfineRed,
                   Eigen::VectorX<ScalarType>& dFineRed) const {
      assemblers[finestLevel].createReducedVector(dFineFull, dFineRed);
      dFineRed = iterativeSolver.solve(RfineRed);
    }

    void transformToFineFull(const Eigen::VectorX<double>& dFineRed,Eigen::VectorX<double>& dFineFull)
    {
      assemblers[finestLevel].createFullVector(dFineRed, dFineFull);
    }

    template <typename ScalarType>
    void solve(Eigen::VectorX<ScalarType>& dFineRed, const Eigen::VectorX<ScalarType>& RfineRed_) const {
      //    Eigen::MatrixXd KcoarseRed,KfineRed;
      Eigen::VectorXd RcoarseRed, dCoarseRed, dCoarseFull,RfineFull,RcoarseFull;
      dCoarseFull.setZero(assemblers[finestLevel-1].size());
      Eigen::VectorXd dFineFull;
      dFineFull.setZero(assemblers[finestLevel].size());

      double lambdaLoad = -1;

      requirementType = Ikarus::FErequirementsBuilder()
                            .insertGlobalSolution(Ikarus::FESolutions::displacement, dFineFull)
                            .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                            .addAffordance(Ikarus::VectorAffordances::forces)
                            .build();
      Eigen::VectorXd RfineRed = -assemblers[finestLevel].getReducedVector(requirementType);
      requirementType   = Ikarus::FErequirementsBuilder()
                            .insertGlobalSolution(Ikarus::FESolutions::displacement, dCoarseFull)
                            .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                            .addAffordance(Ikarus::VectorAffordances::forces)
                            .build();
      auto& KcoarseRed         = assemblers[finestLevel-1].getReducedMatrix(requirementType);

      directSolver.compute(KcoarseRed);
      assemblers[finestLevel].createFullVector(RfineRed, RfineFull);
      transfer.restrictTo(finestLevel-1, RfineFull, RcoarseFull);
      assemblers[finestLevel-1].createReducedVector(RcoarseFull, RcoarseRed);

      directSolver.solve(dCoarseRed, RcoarseRed);

      assemblers[finestLevel-1].createFullVector(dCoarseRed, dCoarseFull);

      transfer.prolongateFrom(finestLevel-1, dCoarseFull, dFineFull);
      requirementType = Ikarus::FErequirementsBuilder()
                            .insertGlobalSolution(Ikarus::FESolutions::displacement, dFineFull)
                            .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                            .addAffordance(Ikarus::VectorAffordances::forces)
                            .build();
      auto& KfineRed = assemblers[finestLevel].getReducedMatrix(requirementType);
      iterativeSolver.compute(KfineRed);
      assemblers[finestLevel].createReducedVector(dFineFull, dFineRed);
      // Pre-Smoothing
      smoothing(dFineFull, RfineRed, dFineRed);

      Eigen::VectorXd residualMGFineRed = RfineRed - KfineRed * dFineRed;
      Eigen::VectorXd residualMGFineFull;
      assemblers[finestLevel].createFullVector(residualMGFineRed, residualMGFineFull);
      Eigen::VectorXd residualMGCoarseFull;
      transfer.restrictTo(finestLevel-1, residualMGFineFull, residualMGCoarseFull);
      Eigen::VectorXd residualMGCoarseRed;
      assemblers[finestLevel-1].createReducedVector(residualMGCoarseFull, residualMGCoarseRed);

      Eigen::VectorXd eCoarseRed;
      directSolver.solve(eCoarseRed, residualMGCoarseRed);
      Eigen::VectorXd eCoarseFull;
      assemblers[finestLevel-1].createFullVector(eCoarseRed, eCoarseFull);
      Eigen::VectorXd eFineFull;
      transfer.prolongateFrom(finestLevel-1, eCoarseFull, eFineFull);

      dFineFull += eFineFull;
      assemblers[finestLevel].createReducedVector(dFineFull, dFineRed);
      // Post smoothing
      smoothing(dFineFull, RfineRed, dFineRed);
    }
  };
}  // namespace Ikarus
