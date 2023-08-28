// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

namespace Ikarus {

  template <typename Basis, typename FEContainer>
  Eigen::VectorXd FlatAssemblerBase<Basis, FEContainer>::createFullVector(Eigen::Ref<const Eigen::VectorXd> reducedVector) {
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
  Eigen::VectorXd &VectorFlatAssembler<Basis, FEContainer>::getRawVectorImpl(const FERequirementType &feRequirements) {
    vecRaw.setZero(this->size());
    Eigen::VectorXd vecLocal;
    std::vector<GlobalIndex> dofs;
    for (auto &fe : this->finiteElements()) {
      vecLocal.setZero(fe.size());
      dofs.resize(0);
      fe.calculateVector(feRequirements, vecLocal);
      fe.globalFlatIndices(dofs);
      for (int i = 0; auto id : dofs) {
        vecRaw(id[0]) += vecLocal(i);
        ++i;
      }
    }
    return vecRaw;
  }

  template <typename Basis, typename FEContainer>
  Eigen::VectorXd &VectorFlatAssembler<Basis, FEContainer>::getVectorImpl(const FERequirementType &feRequirements) {
    vec = getRawVectorImpl(feRequirements);
    for (auto i = 0U; i < this->size(); ++i) {
      if (this->isConstrained(i)) vec[i] = 0;
    }
    return vec;
  }

  template <typename Basis, typename FEContainer>
  Eigen::VectorXd &VectorFlatAssembler<Basis, FEContainer>::getReducedVectorImpl(
      const FERequirementType &feRequirements) {
    vecRed.setZero(this->reducedSize());
    int reducedCounter = 0;
    Eigen::VectorXd vecLocal;
    std::vector<GlobalIndex> dofs;
    for (auto &fe : this->finiteElements()) {
      vecLocal.setZero(fe.size());
      dofs.resize(0);
      fe.calculateVector(feRequirements, vecLocal);
      fe.globalFlatIndices(dofs);
      assert(static_cast<long int>(dofs.size()) == vecLocal.size() && "The returned vector has wrong rowSize!");
      for (int i = 0; auto &&dofIndex : dofs) {
        if (this->isConstrained(dofIndex)) {
          ++reducedCounter;
          ++i;
          continue;
        } else
          vecRed(dofIndex[0] - this->constraintsBelow(dofIndex[0])) += vecLocal[i++];
      }
    }
    return vecRed;
  }

  template <typename Basis, typename FEContainer>
  Eigen::SparseMatrix<double> &SparseFlatAssembler<Basis, FEContainer>::getRawMatrixImpl(
      const FERequirementType &feRequirements) {
    if (!isOccupationPatternCreated) createOccupationPattern();
    if (!areLinearDofsPerElementCreated) createLinearDOFsPerElement();
    spMatRaw.coeffs().setZero();
    Eigen::MatrixXd A;
    for (size_t elementIndex = 0; const auto &fe : this->finiteElements()) {
      A.setZero(fe.size(), fe.size());
      fe.calculateMatrix(feRequirements, A);
      assert(std::sqrt(elementLinearIndices[elementIndex].size()) == A.rows()
             && "The returned matrix has wrong rowSize!");
      assert(std::sqrt(elementLinearIndices[elementIndex].size()) == A.cols()
             && "The returned matrix has wrong colSize!");
      for (Eigen::Index linearIndex = 0; double matrixEntry : A.reshaped())
        spMatRaw.coeffs()(elementLinearIndices[elementIndex][linearIndex++]) += matrixEntry;
      ++elementIndex;
    }
    return spMatRaw;
  }

  template <typename Basis, typename FEContainer>
  Eigen::SparseMatrix<double> &SparseFlatAssembler<Basis, FEContainer>::getMatrixImpl(
      const FERequirementType &feRequirements) {
    spMat = getRawMatrixImpl(feRequirements);
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
    if (!isReducedOccupationPatternCreated) createReducedOccupationPattern();
    if (!areLinearReducedDofsPerElementCreated) createLinearDOFsPerElementReduced();
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
  void SparseFlatAssembler<Basis, FEContainer>::createOccupationPattern() {
    spMatRaw.resize(this->size(), this->size());
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

    spMatRaw.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
    isOccupationPatternCreated = true;
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::createReducedOccupationPattern() {
    spMatReduced.resize(this->reducedSize(), this->reducedSize());
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

    spMatReduced.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
    isReducedOccupationPatternCreated = true;
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::createLinearDOFsPerElement() {
    std::vector<GlobalIndex> dofs;
    for (auto &&fe : this->finiteElements()) {
      dofs.resize(0);
      fe.globalFlatIndices(dofs);
      elementLinearIndices.emplace_back(Dune::power(dofs.size(),2));
      for (Eigen::Index linearIndexOfElement = 0; auto &&c : dofs)
        for (auto &&r : dofs)
          elementLinearIndices.back()[linearIndexOfElement++] = spMatRaw.getLinearIndex(r[0], c[0]);
    }
    areLinearDofsPerElementCreated = true;
  }

  template <typename Basis, typename FEContainer>
  void SparseFlatAssembler<Basis, FEContainer>::createLinearDOFsPerElementReduced() {
    std::vector<GlobalIndex> dofs;
    for (auto &&fe : this->finiteElements()) {
      dofs.resize(0);
      fe.globalFlatIndices(dofs);
      elementLinearReducedIndices.emplace_back();
      for (auto r = 0U; r < dofs.size(); ++r) {
        if (this->isConstrained(dofs[r][0])) continue;
        for (auto c = 0U; c < dofs.size(); ++c) {
          if (this->isConstrained(dofs[c][0])) continue;
          elementLinearReducedIndices.back().push_back(spMatReduced.getLinearIndex(
              dofs[r][0] - this->constraintsBelow(dofs[r][0]), dofs[c][0] - this->constraintsBelow(dofs[c][0])));
        }
      }
    }
    areLinearReducedDofsPerElementCreated = true;
  }

  template <typename Basis, typename FEContainer>
  Eigen::MatrixXd &DenseFlatAssembler<Basis, FEContainer>::getReducedMatrixImpl(const FERequirementType &feRequirements) {
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

  template <typename Basis, typename FEContainer>
  Eigen::MatrixXd &DenseFlatAssembler<Basis, FEContainer>::getRawMatrixImpl(const FERequirementType &feRequirements) {
    matRaw.setZero(this->size(), this->size());
    Eigen::MatrixXd matLocal;
    std::vector<GlobalIndex> dofs;
    for (auto &fe : this->finiteElements()) {
      matLocal.setZero(fe.size(), fe.size());
      dofs.resize(0);
      fe.calculateMatrix(feRequirements, matLocal);
      fe.globalFlatIndices(dofs);
      for (auto i = 0; auto idi : dofs) {
        for (auto j = 0; auto idj : dofs) {
          matRaw(idi[0], idj[0]) += matLocal(i, j);
          ++j;
        }
        ++i;
      }
    }
    return matRaw;
  }

  template <typename Basis, typename FEContainer>
  Eigen::MatrixXd &DenseFlatAssembler<Basis, FEContainer>::getMatrixImpl(const FERequirementType &feRequirements) {
    mat = getRawMatrixImpl(feRequirements);
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) mat.col(i).setZero();
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) mat.row(i).setZero();
    for (auto i = 0U; i < this->size(); ++i)
      if (this->isConstrained(i)) mat(i, i) = 1;
    return mat;
  }

}  // namespace Ikarus
