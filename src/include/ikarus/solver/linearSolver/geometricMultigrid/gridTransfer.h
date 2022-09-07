//
// Created by ac129893 on 21.03.2022.
//

#pragma once
#include <memory>
#include <vector>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Ikarus {
  template <typename Grid>
  class GridTransfer {
    static constexpr int gridDim = Grid::dimension;

  public:
    GridTransfer(const std::shared_ptr<Grid>& p_grid) : grid{p_grid.get()} {}
    GridTransfer(const std::unique_ptr<Grid>& p_grid) : grid{p_grid.get()} {}
    GridTransfer(const Grid* p_grid) : grid{p_grid} {}

    void prolongateFrom(int coarseID, const Eigen::VectorXd& coarse, Eigen::VectorXd& fine) const {
      fine = transferMatrices[coarseID] * coarse;
    }

    void restrictTo(int coarseID, const Eigen::VectorXd& fine, Eigen::VectorXd& coarse) const {
      coarse = transferMatrices[coarseID].transpose() * fine;
    }

    template <typename PreBasisFactory>
    void createOperators(const PreBasisFactory& preBasisFactory) {
      transferMatrices.resize(grid->maxLevel());

      std::vector<Eigen::Triplet<double>> vectorOfTriples;
      //    vectorOfTriples. //Number of Vertices times 8 as estimate
      //    vectorOfTriples.reserve(grid->leafGridView().size(gridDim) * 8); //Number of Vertices times 8 as estimate

      for (int level = 0; level < grid->maxLevel(); ++level) {
        vectorOfTriples.clear();
        auto coarseGridView = grid->levelGridView(level);
        auto fineGridView   = grid->levelGridView(level + 1);

        const auto& coarseIndexSet = coarseGridView.indexSet();
        const auto& fineIndexSet   = fineGridView.indexSet();

        auto coarseBasis            = makeBasis(coarseGridView, preBasisFactory);
        constexpr int numDofPerNode = decltype(coarseBasis)::PreBasis::Node::CHILDREN;
        auto fineBasis              = makeBasis(fineGridView, preBasisFactory);
        auto coarseLocalView        = coarseBasis.localView();
        auto fineLocalView          = fineBasis.localView();

        std::vector<Dune::FieldVector<double, 1>> NcoarseEvaluated;
        std::vector<Dune::FieldVector<double, gridDim>> lagrangeNodeCoords;

        for (auto& coarseElement : elements(coarseGridView)) {
          coarseLocalView.bind(coarseElement);
          const auto& coarseFE = coarseLocalView.tree().child(0).finiteElement();
          const int numNCoarse = coarseFE.localBasis().size();  // Chapter 8<
          NcoarseEvaluated.resize(numNCoarse);

          for (auto& childsElement : descendantElements(coarseElement, coarseElement.level() + 1)) {
            fineLocalView.bind(childsElement);
            const auto& fineFE = fineLocalView.tree().child(0).finiteElement();
            const int numNFine = fineFE.localBasis().size();

            obtainLagrangeNodePositions(fineFE, lagrangeNodeCoords);  // CoarseIndex Set Chapter 5.6

            const auto geoInFather = childsElement.geometryInFather();
            for (int i = 0; i < numNFine; ++i) {
              const auto localInFather = geoInFather.global(lagrangeNodeCoords[i]);
              coarseFE.localBasis().evaluateFunction(localInFather, NcoarseEvaluated);

              for (int j = 0; j < numNCoarse; ++j) {
                for (int k = 0; k < numDofPerNode; ++k) {
                  const size_t globalFine   = fineLocalView.index((fineLocalView.tree().child(k).localIndex(i)));
                  const size_t globalCoarse = coarseLocalView.index((coarseLocalView.tree().child(k).localIndex(j)));
                  vectorOfTriples.emplace_back(globalFine, globalCoarse, NcoarseEvaluated[j]);
                }
              }
            }
          }
        }
        transferMatrices[level].resize(fineBasis.size(), coarseBasis.size());
        // create sparse matrix and ignore duplicates
        transferMatrices[level].setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end(),
                                                [](const double&, const double& b) { return b; });
      }
    }

  private:
    template <typename LocalFE>  // Dune Book Page 314
    void obtainLagrangeNodePositions(const LocalFE& localFE,
                                     std::vector<Dune::FieldVector<double, gridDim>>& lagrangeNodeCoords) {
      lagrangeNodeCoords.resize(localFE.size());
      std::vector<double> out;
      for (int i = 0; i < gridDim; i++) {
        auto ithCoord = [&i](const Dune::FieldVector<double, gridDim>& x) { return x[i]; };

        localFE.localInterpolation().interpolate(ithCoord, out);

        for (std::size_t j = 0; j < out.size(); j++)
          lagrangeNodeCoords[j][i] = out[j];
      }
    }

    std::vector<Eigen::SparseMatrix<double>> transferMatrices;
    const Grid* grid;
  };

}  // namespace Ikarus
