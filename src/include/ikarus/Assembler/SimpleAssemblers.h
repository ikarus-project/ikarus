//
// Created by Alex on 26.06.2021.
//

#pragma once
#include <ikarus/FiniteElements/FiniteElementPolicies.h>

#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Ikarus::Assembler {
  template <typename DofManagerType>
  class VectorAssembler {
  public:
    explicit VectorAssembler(DofManagerType& dofManager)
        : dofManager_{&dofManager}, vec{Eigen::VectorXd::Zero(dofManager_->correctionSize())} {}

    Eigen::VectorXd& getVector(Ikarus::FiniteElements::VectorAffordances VectorAffordances) {
      return getVectorImpl(VectorAffordances);
    }

  private:
    Eigen::VectorXd& getVectorImpl(Ikarus::FiniteElements::VectorAffordances vecAffordances) {
      vec.setZero(dofManager_->correctionSize());
      for (auto [fe, dofs, vars] : dofManager_->elementDofsVariableTuple()) {
        assert(dofs.size() == calculateVector(fe, vars, vecAffordances).size()
               && "The returned vector has wrong rowSize!");
        vec(dofs) += calculateVector(fe,vars, vecAffordances);
      }
      return vec;
    }

  private:
    DofManagerType* dofManager_;
    Eigen::VectorXd vec{};
  };

  template <typename DofManagerType>
  class DenseMatrixAssembler {
  public:
    explicit DenseMatrixAssembler(DofManagerType& dofManager)
        : dofManager_{&dofManager},
          mat{Eigen::MatrixXd::Zero(dofManager_->correctionSize(), dofManager_->correctionSize())} {}

    Eigen::MatrixXd& getMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances) {
      return getMatrixImpl(MatrixAffordances);
    }

  private:
    Eigen::MatrixXd& getMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matAffordances) {
      mat.setZero(dofManager_->correctionSize(), dofManager_->correctionSize());
      for (auto [fe, dofs, vars] : dofManager_->elementDofsVariableTuple()) {
        assert(dofs.size() == calculateMatrix(fe, vars, matAffordances).rows()
               && "The returned matrix has wrong rowSize!");
        assert(dofs.size() == calculateMatrix(fe, vars, matAffordances).cols()
               && "The returned matrix has wrong colSize!");
        mat(dofs, dofs) += calculateMatrix(fe, vars, matAffordances);
      }
      return mat;
    }

  private:
    DofManagerType* dofManager_;
    Eigen::MatrixXd mat{};
  };

  template <typename DofManagerType>
  class SparseMatrixAssembler {
  public:
    explicit SparseMatrixAssembler(DofManagerType& dofManager) : dofManager_{&dofManager} {}

    Eigen::SparseMatrix<double>& getMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances) {
      return getMatrixImpl(MatrixAffordances);
    }

  private:
    Eigen::SparseMatrix<double>& getMatrixImpl(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances) {
      if (!isOccupationPatternCreated) createOccupationPattern();
      if (!arelinearDofsPerElementCreated) createlinearDofsPerElement();
      spMat.coeffs().setZero();
      Eigen::MatrixXd A;
      for (size_t elementIndex = 0; auto [fe, dofs, vars] : dofManager_->elementDofsVariableTuple()) {
        A = calculateMatrix(fe, vars, MatrixAffordances);
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
      spMat.resize(dofManager_->correctionSize(), dofManager_->correctionSize());
      std::vector<Eigen::Triplet<double>> vectorOfTriples;
      int estimateOfConnectivity = 6;
      vectorOfTriples.reserve(estimateOfConnectivity * vertices((*dofManager_->getGridView())).size());
      for (auto&& dofsOfElement : dofManager_->elementDofs())
        for (auto&& c : dofsOfElement)
          for (auto&& r : dofsOfElement)
            vectorOfTriples.emplace_back(r, c, 0.0);

      spMat.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
      isOccupationPatternCreated = true;
    }

    // This function save the indices of each element in the underlying vector which stores the sparse matrix entries
    void createlinearDofsPerElement() {
      for (Eigen::Index eleIndex = 0; auto&& dofsOfElement : dofManager_->elementDofs()) {
        elementLinearIndices.emplace_back(dofsOfElement.size() * dofsOfElement.size());
        for (Eigen::Index linearIndexOfElement = 0; auto&& c : dofsOfElement)
          for (auto&& r : dofsOfElement)
            elementLinearIndices[eleIndex](linearIndexOfElement++) = spMat.getLinearIndex(r, c);
        ++eleIndex;
      }
      arelinearDofsPerElementCreated = true;
    }

  private:
    bool isOccupationPatternCreated{false};
    bool arelinearDofsPerElementCreated{false};
    DofManagerType* dofManager_;
    Eigen::SparseMatrix<double> spMat;
    std::vector<Eigen::ArrayX<Eigen::Index>> elementLinearIndices;
  };
}  // namespace Ikarus::Assembler