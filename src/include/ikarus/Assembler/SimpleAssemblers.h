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

  template <typename DofManagerType>
  class ScalarAssembler {
  public:
    explicit ScalarAssembler(DofManagerType& dofManager) : feManager_{&dofManager} {}

    double& getScalar(Ikarus::FiniteElements::ScalarAffordances scalarAffordances) {
      return getScalarImpl(scalarAffordances);
    }

  private:
    double& getScalarImpl(Ikarus::FiniteElements::ScalarAffordances scalarAffordances) {
      scal = 0.0;
      for (auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements feParameter;
        feParameter.scalarAffordances = scalarAffordances;
        feParameter.variables         = vars;
        scal += calculateScalar(fe, feParameter);
      }
      return scal;
    }

    DofManagerType* feManager_;
    double scal{0.0};
  };

  template <typename DofManagerType>
  class VectorAssembler {
  public:
    explicit VectorAssembler(DofManagerType& dofManager)
        : feManager_{&dofManager}, vec{Eigen::VectorXd::Zero(feManager_->numberOfDegreesOfFreedom())} {}

    Eigen::VectorXd& getVector(Ikarus::FiniteElements::VectorAffordances vectorAffordances) {
      return getVectorImpl(vectorAffordances);
    }

  private:
    Eigen::VectorXd& getVectorImpl(Ikarus::FiniteElements::VectorAffordances vecAffordances) {
      vec.setZero(feManager_->numberOfDegreesOfFreedom());
      for (auto&& [fe, dofIndices, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements feParameter;
        feParameter.vectorAffordances = vecAffordances;
        feParameter.variables         = vars;
        assert(dofIndices.size() == calculateVector(fe, feParameter).size()
               && "The returned vector has wrong rowSize!");
        vec(dofIndices) += calculateVector(fe, feParameter);
      }
      return vec;
    }

    DofManagerType* feManager_;
    Eigen::VectorXd vec{};
  };

  template <typename DofManagerType>
  class DenseMatrixAssembler {
  public:
    explicit DenseMatrixAssembler(DofManagerType& dofManager)
        : feManager_{&dofManager},
          mat{Eigen::MatrixXd::Zero(feManager_->numberOfDegreesOfFreedom(), feManager_->numberOfDegreesOfFreedom())} {}

    Eigen::MatrixXd& getMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances) {
      return getMatrixImpl(MatrixAffordances);
    }

  private:
    Eigen::MatrixXd& getMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matAffordances) {
      mat.setZero(feManager_->numberOfDegreesOfFreedom(), feManager_->numberOfDegreesOfFreedom());
      for (auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements feParameter;
        feParameter.matrixAffordances = matAffordances;
        feParameter.variables         = vars;
        assert(dofs.size() == calculateMatrix(fe, feParameter).rows()
               && "The returned matrix has wrong rowSize!");
        assert(dofs.size() == calculateMatrix(fe, feParameter).cols()
               && "The returned matrix has wrong colSize!");
        mat(dofs, dofs) += calculateMatrix(fe, feParameter);
      }
      return mat;
    }

    DofManagerType* feManager_;
    Eigen::MatrixXd mat{};
  };

  template <typename DofManagerType>
  class SparseMatrixAssembler {
  public:
    explicit SparseMatrixAssembler(DofManagerType& dofManager) : feManager_{&dofManager} {}

    Eigen::SparseMatrix<double>& getMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances) {
      return getMatrixImpl(MatrixAffordances);
    }

  private:
    Eigen::SparseMatrix<double>& getMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matrixAffordances) {
      if (!isOccupationPatternCreated) createOccupationPattern();
      if (!arelinearDofsPerElementCreated) createlinearDofsPerElement();
      spMat.coeffs().setZero();
      Eigen::MatrixXd A;
      for (size_t elementIndex = 0; auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
        FiniteElements::FErequirements feParameter;
        feParameter.matrixAffordances = matrixAffordances;
        feParameter.variables         = vars;
        A = calculateMatrix(fe, feParameter);
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
    DofManagerType* feManager_;
    Eigen::SparseMatrix<double> spMat;
    std::vector<Eigen::ArrayX<Eigen::Index>> elementLinearIndices;
  };
}  // namespace Ikarus::Assembler