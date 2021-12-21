//
// Created by Alex on 26.06.2021.
//

#pragma once
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>

#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"
#include <utility>

#include <dune/common/power.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Ikarus::Assembler {

  template <typename FEManager>
  class ScalarAssembler {
  public:
    explicit ScalarAssembler(FEManager& dofManager) : feManager_{&dofManager} {}

    double& getScalar(Ikarus::FiniteElements::ScalarAffordances scalarAffordances,
                      const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      return getScalarImpl(scalarAffordances, feParameter);
    }

  private:
    double& getScalarImpl(Ikarus::FiniteElements::ScalarAffordances scalarAffordances,
                          const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      scal = 0.0;
      for (auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements fErequirements;
        fErequirements.scalarAffordances = scalarAffordances;
        fErequirements.variables         = vars;
        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
        scal += calculateScalar(fe, fErequirements);
      }
      return scal;
    }

    FEManager* feManager_;
    double scal{0.0};
  };

  template <typename FEManager, typename DirichletManager>
  class VectorAssembler {
  public:
    explicit VectorAssembler(const FEManager& dofManager, const DirichletManager& dirichletManager)
        : feManager_{&dofManager}, dirichletManager_{&dirichletManager} {}

    Eigen::VectorXd& getVector(Ikarus::FiniteElements::VectorAffordances vectorAffordances,
                               const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      return getVectorImpl(vectorAffordances, feParameter);
    }

    Eigen::VectorXd& getReducedVector(Ikarus::FiniteElements::VectorAffordances vectorAffordances,
                                      const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      return getReducedVectorImpl(vectorAffordances, feParameter);
    }

    auto createFullVector(const Eigen::VectorXd& reducedVector) {
      assert(reducedVector.size() == static_cast<Eigen::Index>(dirichletManager_->numberOfReducedDegreesOfFreedom())
             && "The reduced vector you passed has the wrong dimensions.");
      Eigen::Index reducedCounter = 0;
      Eigen::VectorXd fullVec(feManager_->numberOfDegreesOfFreedom());
      for (Eigen::Index i = 0; i < fullVec.size(); ++i) {
        if (dirichletManager_->isConstrained(i)) {
          ++reducedCounter;
          fullVec[i] = 0.0;
          continue;
        } else
          fullVec[i] = reducedVector[i - reducedCounter];
      }
      return fullVec;
    }

  private:
    Eigen::VectorXd& getVectorImpl(Ikarus::FiniteElements::VectorAffordances vecAffordances,
                                   const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      vec.setZero(feManager_->numberOfDegreesOfFreedom());
      for (auto&& [fe, dofIndices, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements fErequirements;
        fErequirements.vectorAffordances = vecAffordances;
        fErequirements.variables         = vars;
        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
        assert(dofIndices.size() == calculateVector(fe, fErequirements).size()
               && "The returned vector has wrong rowSize!");
        vec(dofIndices) += calculateVector(fe, fErequirements);
      }
      return vec;
    }

    Eigen::VectorXd& getReducedVectorImpl(Ikarus::FiniteElements::VectorAffordances vecAffordances,
                                          const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      vecRed.setZero(dirichletManager_->numberOfReducedDegreesOfFreedom());
      int reducedCounter = 0;
      for (auto&& [fe, dofIndices, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements fErequirements;
        fErequirements.vectorAffordances = vecAffordances;
        fErequirements.variables         = vars;
        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});

        const auto f = calculateVector(fe, fErequirements);
        assert(dofIndices.size() == f.size() && "The returned vector has wrong rowSize!");
        for (int i = 0; auto&& dofIndex : dofIndices) {
          if (dirichletManager_->isConstrained(dofIndex)) {
            ++reducedCounter;
            ++i;
            continue;
          } else
            vecRed(dofIndex - reducedCounter) += f[i++];
        }
      }
      return vecRed;
    }

    FEManager const* feManager_;
    DirichletManager const* dirichletManager_;
    Eigen::VectorXd vec{};
    Eigen::VectorXd vecRed{};
  };

  template <typename FEManager, typename DirichletManager>
  class DenseMatrixAssembler {
  public:
    explicit DenseMatrixAssembler(FEManager& dofManager, const DirichletManager& dirichletManager)
        : feManager_{&dofManager}, dirichletManager_{&dirichletManager} {}

    Eigen::MatrixXd& getMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances,
                               const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      return getMatrixImpl(MatrixAffordances, feParameter);
    }

    Eigen::MatrixXd& getReducedMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances,
                                      const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      return getReducedMatrixImpl(MatrixAffordances, feParameter);
    }

  private:
    Eigen::MatrixXd& getMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matAffordances,
                                   const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      mat.setZero(feManager_->numberOfDegreesOfFreedom(), feManager_->numberOfDegreesOfFreedom());
      for (auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements fErequirements;
        fErequirements.matrixAffordances = matAffordances;
        fErequirements.variables         = vars;
        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
        assert(dofs.size() == calculateMatrix(fe, fErequirements).rows() && "The returned matrix has wrong rowSize!");
        assert(dofs.size() == calculateMatrix(fe, fErequirements).cols() && "The returned matrix has wrong colSize!");
        mat(dofs, dofs) += calculateMatrix(fe, fErequirements);
      }
      return mat;
    }

    Eigen::MatrixXd& getReducedMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matAffordances,
                                          const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      matRed.setZero(dirichletManager_->numberOfReducedDegreesOfFreedom(),
                     dirichletManager_->numberOfReducedDegreesOfFreedom());
      for (auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements fErequirements;
        fErequirements.matrixAffordances = matAffordances;
        fErequirements.variables         = vars;
        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
        const auto eleMat = calculateMatrix(fe, fErequirements);
        assert(dofs.size() == eleMat.rows() && "The returned matrix has wrong rowSize!");
        assert(dofs.size() == eleMat.cols() && "The returned matrix has wrong colSize!");
        for (int r = 0; r < dofs.size(); ++r) {
          if (dirichletManager_->isConstrained(dofs[r]))
            continue;
          else {
            for (int c = 0; c < dofs.size(); ++c) {
              if (dirichletManager_->isConstrained(dofs[c])) continue;
              matRed(dofs[r]-dirichletManager_->constraintsBelow(dofs[r]), dofs[c]-dirichletManager_->constraintsBelow(dofs[c])) += eleMat(r, c);
            }
          }
        }
      }
      return matRed;
    }

    FEManager* feManager_;
    DirichletManager const* dirichletManager_;
    Eigen::MatrixXd mat{};
    Eigen::MatrixXd matRed{};
  };

  template <typename FEManager>
  class SparseMatrixAssembler {
  public:
    explicit SparseMatrixAssembler(FEManager& dofManager) : feManager_{&dofManager} {}

    Eigen::SparseMatrix<double>& getMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances,
                                           const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      return getMatrixImpl(MatrixAffordances, feParameter);
    }

    Eigen::SparseMatrix<double>& getMatrixReduced(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances,
                                                  const std::optional<FEParameterValuePair>& feParameter
                                                  = std::nullopt) {
      return getMatrixImpl(MatrixAffordances, feParameter);
    }

  private:
    Eigen::SparseMatrix<double>& getMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matrixAffordances,
                                               const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
      if (!isOccupationPatternCreated) createOccupationPattern();
      if (!arelinearDofsPerElementCreated) createlinearDofsPerElement();
      spMat.coeffs().setZero();
      Eigen::MatrixXd A;
      for (size_t elementIndex = 0; auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements fErequirements;
        fErequirements.matrixAffordances = matrixAffordances;
        fErequirements.variables         = vars;
        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
        A = calculateMatrix(fe, fErequirements);
        assert(dofs.size() == A.rows() && "The returned matrix has wrong rowSize!");
        assert(dofs.size() == A.cols() && "The returned matrix has wrong colSize!");
        for (Eigen::Index linearIndex = 0; double matrixEntry : A.reshaped())
          spMat.coeffs()(elementLinearIndices[elementIndex](linearIndex++)) += matrixEntry;
        ++elementIndex;
      }
      return spMat;
    }

    // https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
    void createOccupationPattern() {
      spMat.resize(feManager_->numberOfDegreesOfFreedom(), feManager_->numberOfDegreesOfFreedom());
      std::vector<Eigen::Triplet<double>> vectorOfTriples;
      int estimateOfConnectivity = 6;
      vectorOfTriples.reserve(estimateOfConnectivity * vertices((*feManager_->getGridView())).size());
      for (auto&& dofsOfElement : feManager_->elementDofs())
        for (auto&& c : dofsOfElement)
          for (auto&& r : dofsOfElement)
            vectorOfTriples.emplace_back(r, c, 0.0);

      spMat.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
      isOccupationPatternCreated = true;
    }

    // This function save the indices of each element in the underlying vector which stores the sparse matrix entries
    void createlinearDofsPerElement() {
      for (Eigen::Index eleIndex = 0; auto&& dofsOfElement : feManager_->elementDofs()) {
        elementLinearIndices.emplace_back(Dune::Power<2>::eval(dofsOfElement.size()));
        for (Eigen::Index linearIndexOfElement = 0; auto&& c : dofsOfElement)
          for (auto&& r : dofsOfElement)
            elementLinearIndices[eleIndex](linearIndexOfElement++) = spMat.getLinearIndex(r, c);
        ++eleIndex;
      }
    }

    bool isOccupationPatternCreated{false};
    bool arelinearDofsPerElementCreated{false};
    FEManager* feManager_;
    Eigen::SparseMatrix<double> spMat;
    std::vector<Eigen::ArrayX<Eigen::Index>> elementLinearIndices;
  };
}  // namespace Ikarus::Assembler