// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

namespace Ikarus {

  template <typename Basis, typename FEContainer>
  Eigen::VectorXd FlatAssemblerBase<Basis, FEContainer>::createFullVector(
      Eigen::Ref<const Eigen::VectorXd> reducedVector) {
    assert(reducedVector.size() == static_cast<Eigen::Index>(this->reducedSize())
           && "The reduced vector you passed has the wrong dimensions.");
    Eigen::Index reducedCounter = 0;
    Eigen::VectorXd fullVec(this->size());
    for (Eigen::Index i = 0; i < fullVec.size(); ++i) {
      if (dirichletValues->isConstrained(i)) {
        ++reducedCounter;
        fullVec[i] = 0.0;
        continue;
      } else
        fullVec[i] = reducedVector[i - reducedCounter];
    }
    return fullVec;
  }

  template <typename Basis, typename FEContainer>
  void VectorFlatAssembler<Basis, FEContainer>::assembleRawVectorImpl(const FERequirementType &feRequirements,
                                                                      Eigen::VectorXd &assemblyVec) {
    assemblyVec.setZero(this->size());
    Eigen::VectorXd vecLocal;
    std::vector<GlobalIndex> dofs;
    for (auto &fe : this->finiteElements()) {
      vecLocal.setZero(fe.size());
      dofs.resize(0);
      fe.calculateVector(feRequirements, vecLocal);
      fe.globalFlatIndices(dofs);
      for (int i = 0; auto id : dofs) {
        assemblyVec(id[0]) += vecLocal(i);
        ++i;
      }
    }
  }

  template <typename Basis, typename FEContainer>
  Eigen::VectorXd &VectorFlatAssembler<Basis, FEContainer>::getRawVectorImpl(const FERequirementType &feRequirements) {
    assembleRawVectorImpl(feRequirements, vecRaw);
    return vecRaw;
  }

  template <typename Basis, typename FEContainer>
  Eigen::VectorXd &VectorFlatAssembler<Basis, FEContainer>::getVectorImpl(const FERequirementType &feRequirements) {
    assembleRawVectorImpl(feRequirements, vec);
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) vec[i] = 0;
    return vec;
  }

  template <typename Basis, typename FEContainer>
  Eigen::VectorXd &VectorFlatAssembler<Basis, FEContainer>::getReducedVectorImpl(
      const FERequirementType &feRequirements) {
    vecRed.setZero(this->reducedSize());
    Eigen::VectorXd vecLocal;
    std::vector<GlobalIndex> dofs;
    for ( auto &fe : this->finiteElements()) {
      vecLocal.setZero(fe.size());
      dofs.resize(0);
      fe.calculateVector(feRequirements, vecLocal);
      fe.globalFlatIndices(dofs);
      assert(static_cast<long int>(dofs.size()) == vecLocal.size() && "The returned vector has wrong rowSize!");
      for (int localConstrainedCounter = 0; auto &&dofIndex : dofs) {
        if (this->isConstrained(dofIndex)) {
          ++localConstrainedCounter;
          continue;
        } else
          vecRed(dofIndex[0] - this->constraintsBelow(dofIndex[0])) += vecLocal[localConstrainedCounter++];
      }
    }
    return vecRed;
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::assembleRawMatrixImpl(const FERequirementType &feRequirements,
                                                                      Eigen::SparseMatrix<double> &assemblyMat) {
    assemblyMat.coeffs().setZero();
    Eigen::MatrixXd A;
    for (size_t elementIndex = 0; const auto &fe : this->finiteElements()) {
      A.setZero(fe.size(), fe.size());
      fe.calculateMatrix(feRequirements, A);
      assert(static_cast<Eigen::Index>(std::sqrt(elementLinearIndices[elementIndex].size())) == A.rows()
             && "The returned matrix has wrong rowSize!");
      assert(static_cast<Eigen::Index>(std::sqrt(elementLinearIndices[elementIndex].size())) == A.cols()
             && "The returned matrix has wrong colSize!");
      for (Eigen::Index linearIndex = 0; double matrixEntry : A.reshaped())
        assemblyMat.coeffs()(elementLinearIndices[elementIndex][linearIndex++]) += matrixEntry;
      ++elementIndex;
    }
  }

  template <typename Basis, typename FEContainer>
  Eigen::SparseMatrix<double> &SparseFlatAssembler<Basis, FEContainer>::getRawMatrixImpl(
      const FERequirementType &feRequirements) {
    std::call_once(sparsePreProcessorRaw, [&]() { preProcessSparseMatrix(spMatRaw); });
    assembleRawMatrixImpl(feRequirements, spMatRaw);
    return spMatRaw;
  }

  template <typename Basis, typename FEContainer>
  Eigen::SparseMatrix<double> &SparseFlatAssembler<Basis, FEContainer>::getMatrixImpl(
      const FERequirementType &feRequirements) {
    std::call_once(sparsePreProcessor, [&]() { preProcessSparseMatrix(spMat); });
    assembleRawMatrixImpl(feRequirements, spMat);
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) spMat.col(i) *= 0;
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) spMat.row(i) *= 0;
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) spMat.diagonal()[i] = 1;
    return spMat;
  }

  template <typename Basis, typename FEContainer>
  Eigen::SparseMatrix<double> &SparseFlatAssembler<Basis, FEContainer>::getReducedMatrixImpl(
      const FERequirementType &feRequirements) {
    std::call_once(sparsePreProcessorReduced, [&]() { preProcessSparseMatrixReduced(spMatReduced); });
    spMatReduced.coeffs().setZero();
    Eigen::MatrixXd A;
    std::vector<GlobalIndex> dofs;
    for (size_t elementIndex = 0; const auto &fe : this->finiteElements()) {
      A.setZero(fe.size(), fe.size());
      dofs.resize(0);
      fe.calculateMatrix(feRequirements, A);
      fe.globalFlatIndices(dofs);
      assert(dofs.size() == static_cast<unsigned>(A.rows()) && "The returned matrix has wrong rowSize!");
      assert(dofs.size() == static_cast<unsigned>(A.cols()) && "The returned matrix has wrong colSize!");
      Eigen::Index linearIndex = 0;
      for (auto r = 0U; r < dofs.size(); ++r) {
        if (this->isConstrained(dofs[r][0]))
          continue;
        else {
          for (auto c = 0U; c < dofs.size(); ++c) {
            if (this->isConstrained(dofs[c][0])) continue;
            spMatReduced.coeffs()(elementLinearReducedIndices[elementIndex][linearIndex++]) += A(r, c);
          }
        }
      }
      ++elementIndex;
    }
    return spMatReduced;
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::createOccupationPattern(Eigen::SparseMatrix<double> &assemblyMat) {
    assemblyMat.resize(this->size(), this->size());
    std::vector<Eigen::Triplet<double>> vectorOfTriples;

    vectorOfTriples.reserve(this->estimateOfConnectivity());
    std::vector<GlobalIndex> dofs;
    for (auto &&fe : this->finiteElements()) {
      dofs.resize(0);
      fe.globalFlatIndices(dofs);
      for (auto idi : dofs)
        for (auto idj : dofs)
          vectorOfTriples.emplace_back(idi[0], idj[0], 0.0);
    }
    assemblyMat.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::createReducedOccupationPattern(
      Eigen::SparseMatrix<double> &assemblyMat) {
    assemblyMat.resize(this->reducedSize(), this->reducedSize());
    std::vector<Eigen::Triplet<double>> vectorOfTriples;
    using std::size;

    vectorOfTriples.reserve(this->estimateOfConnectivity());
    std::vector<GlobalIndex> dofs;
    for (auto &fe : this->finiteElements()) {
      dofs.resize(0);
      fe.globalFlatIndices(dofs);
      for (auto r = 0U; r < dofs.size(); ++r) {
        if (this->isConstrained(dofs[r][0]))
          continue;
        else {
          for (auto c = 0U; c < dofs.size(); ++c) {
            if (this->isConstrained(dofs[c][0])) continue;
            vectorOfTriples.emplace_back(dofs[r][0] - this->constraintsBelow(dofs[r][0]),
                                         dofs[c][0] - this->constraintsBelow(dofs[c][0]), 0.0);
          }
        }
      }
    }
    assemblyMat.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::createLinearDOFsPerElement(Eigen::SparseMatrix<double> &assemblyMat) {
    std::vector<GlobalIndex> dofs;
    for (auto &&fe : this->finiteElements()) {
      dofs.resize(0);
      fe.globalFlatIndices(dofs);
      elementLinearIndices.emplace_back(Dune::power(dofs.size(), 2));
      for (Eigen::Index linearIndexOfElement = 0; auto &&c : dofs)
        for (auto &&r : dofs)
          elementLinearIndices.back()[linearIndexOfElement++] = assemblyMat.getLinearIndex(r[0], c[0]);
    }
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::createLinearDOFsPerElementReduced(
      Eigen::SparseMatrix<double> &assemblyMat) {
    std::vector<GlobalIndex> dofs;
    for (auto &&fe : this->finiteElements()) {
      dofs.resize(0);
      fe.globalFlatIndices(dofs);
      elementLinearReducedIndices.emplace_back();
      for (auto r = 0U; r < dofs.size(); ++r) {
        if (this->isConstrained(dofs[r][0])) continue;
        for (auto c = 0U; c < dofs.size(); ++c) {
          if (this->isConstrained(dofs[c][0])) continue;
          elementLinearReducedIndices.back().push_back(assemblyMat.getLinearIndex(
              dofs[r][0] - this->constraintsBelow(dofs[r][0]), dofs[c][0] - this->constraintsBelow(dofs[c][0])));
        }
      }
    }
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::preProcessSparseMatrix(Eigen::SparseMatrix<double> &assemblyMat) {
    createOccupationPattern(assemblyMat);
    createLinearDOFsPerElement(assemblyMat);
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::preProcessSparseMatrixReduced(Eigen::SparseMatrix<double> &assemblyMat) {
    createReducedOccupationPattern(assemblyMat);
    createLinearDOFsPerElementReduced(assemblyMat);
  }

  template <typename Basis, typename FEContainer>
  void DenseFlatAssembler<Basis, FEContainer>::assembleRawMatrixImpl(const FERequirementType &feRequirements,
                                                                     Eigen::MatrixXd &assemblyMat) {
    assemblyMat.setZero(this->size(), this->size());
    Eigen::MatrixXd matLocal;
    std::vector<GlobalIndex> dofs;
    for (auto &fe : this->finiteElements()) {
      matLocal.setZero(fe.size(), fe.size());
      dofs.resize(0);
      fe.calculateMatrix(feRequirements, matLocal);
      fe.globalFlatIndices(dofs);
      for (auto i = 0; auto idi : dofs) {
        for (auto j = 0; auto idj : dofs) {
          assemblyMat(idi[0], idj[0]) += matLocal(i, j);
          ++j;
        }
        ++i;
      }
    }
  }

  template <typename Basis, typename FEContainer>
  Eigen::MatrixXd &DenseFlatAssembler<Basis, FEContainer>::getRawMatrixImpl(const FERequirementType &feRequirements) {
    assembleRawMatrixImpl(feRequirements, matRaw);
    return matRaw;
  }

  template <typename Basis, typename FEContainer>
  Eigen::MatrixXd &DenseFlatAssembler<Basis, FEContainer>::getMatrixImpl(const FERequirementType &feRequirements) {
    assembleRawMatrixImpl(feRequirements, mat);
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) mat.col(i).setZero();
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) mat.row(i).setZero();
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) mat(i, i) = 1;
    return mat;
  }

  template <typename Basis, typename FEContainer>
  Eigen::MatrixXd &DenseFlatAssembler<Basis, FEContainer>::getReducedMatrixImpl(
      const FERequirementType &feRequirements) {
    matRed.setZero(this->reducedSize(), this->reducedSize());
    Eigen::MatrixXd matLocal;
    std::vector<GlobalIndex> dofs;
    for (auto &fe : this->finiteElements()) {
      matLocal.setZero(fe.size(), fe.size());
      dofs.resize(0);
      fe.calculateMatrix(feRequirements, matLocal);
      fe.globalFlatIndices(dofs);
      assert(dofs.size() == static_cast<unsigned>(matLocal.rows()) && "The returned matrix has wrong rowSize!");
      assert(dofs.size() == static_cast<unsigned>(matLocal.cols()) && "The returned matrix has wrong colSize!");
      for (auto r = 0U; r < dofs.size(); ++r) {
        if (this->isConstrained(dofs[r][0])) {
          continue;
        } else {
          for (auto c = 0U; c < dofs.size(); ++c) {
            if (this->isConstrained(dofs[c][0])) {
              continue;
            }
            matRed(dofs[r][0] - this->constraintsBelow(dofs[r][0]), dofs[c][0] - this->constraintsBelow(dofs[c][0]))
                += matLocal(r, c);
          }
        }
      }
    }
    return matRed;
  }
}  // namespace Ikarus
