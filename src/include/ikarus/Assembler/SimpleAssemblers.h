//
// Created by Alex on 26.06.2021.
//

#pragma once
#include "ikarus/FiniteElements/Interface/FiniteElementFunctionConcepts.h"
#include <ikarus/utils/concepts.h>

#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"
#include <ranges>
#include <utility>

#include <dune/common/power.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Ikarus {

  /*!
   * The FlatAssemblerBase takes care of common subtask done by flat assemblers
   */
  template <typename Basis, typename FEContainer>
  class FlatAssemblerBase {
  public:
    using GridView = typename Basis::GridView;
    FlatAssemblerBase(const Basis &basis, const FEContainer &fes, const std::vector<bool> &dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {
      constraintsBelow_.reserve(basis_->size());
      size_t counter = 0;
      for (auto iv : std::ranges::iota_view{size_t(0), basis_->size()}) {
        constraintsBelow_.emplace_back(counter);
        if (dirichFlags[iv]) ++counter;
      }
      fixedDofs = std::ranges::count(dirichFlags, true);
    }

    /**  Returns the size of the free degrees of freeedom, which are not fixed by a dirichlet boundary condition */
    size_t reducedSize() { return basis_->size() - fixedDofs; }

    /**  Returns the size of nodes, i.e. the number of degrees of freedom */
    size_t size() { return basis_->size(); }

    /**  Creates a the fullsized vector of size #Dof and inserts the values of reduced Vector at the "free" degrees
     * of freedom and
     * writes a zero for the fixed doffs */
    auto createFullVector(const Eigen::VectorXd &reducedVector) {
      assert(reducedVector.size() == static_cast<Eigen::Index>(this->reducedSize())
             && "The reduced vector you passed has the wrong dimensions.");
      Eigen::Index reducedCounter = 0;
      Eigen::VectorXd fullVec(this->size());
      for (Eigen::Index i = 0; i < fullVec.size(); ++i) {
        if (dirichletFlags->at(i)) {
          ++reducedCounter;
          fullVec[i] = 0.0;
          continue;
        } else
          fullVec[i] = reducedVector[i - reducedCounter];
      }
      return fullVec;
    }

    /**  Returns the container of finite elements */
    auto &finiteElements() const { return feContainer; }

    /**  Returns the number of constraints below a given degrees of freedom index */
    size_t constraintsBelow(size_t i) const { return constraintsBelow_[i]; }

    /**  Returns the boolean if a given degree of freedom is fixed */
    bool isConstrained(size_t i) const { return dirichletFlags->at(i); }

    /**  Coarse estimate of node cnnectivit, i.e. this relates the bandwith of an sparse matrix.
     * This estimate is désigned that it overestimates the real connectivity since it should
     * only be used fore allocating vectors */
    size_t estimateOfConnectivity() const { return basis_->gridView().size(GridView::dimension) * 8; }

  private:
    Basis const *basis_;
    FEContainer const &feContainer;
    std::vector<bool> const *dirichletFlags;
    std::vector<size_t> constraintsBelow_{};
    size_t fixedDofs{};
  };

  /** ScalarAssembler assembles scalar quantities */
  template <typename Basis, typename FEContainer>
  class ScalarAssembler : public FlatAssemblerBase<Basis, FEContainer> {
    using RequirementType = typename FEContainer::value_type::FERequirementType;

  public:
    ScalarAssembler(const Basis &basis, const FEContainer &fes, const std::vector<bool> &dirichFlags)
        : FlatAssemblerBase<Basis, FEContainer>(basis, fes, dirichFlags) {}

    /** Calculates the scalar quantity which is requested by fErequirements and returns a reference */
    double &getScalar(const RequirementType &fErequirements) { return getScalarImpl(fErequirements); }

  private:
    double &getScalarImpl(const RequirementType &fErequirements) {
      scal = 0.0;
      for (auto &fe : this->finiteElements())
        scal += fe.calculateScalar(fErequirements);
      return scal;
    }

    double scal{0.0};
  };

  /** VectorFlatAssembler assembles vector quantities using a flat basis Indexing strategy */
  template <typename Basis, typename FEContainer>
  class VectorFlatAssembler : public ScalarAssembler<Basis, FEContainer> {
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    using GlobalIndex     = typename FEContainer::value_type::GlobalIndex;

  public:
    VectorFlatAssembler(const Basis &basis, const FEContainer &fes, const std::vector<bool> &dirichFlags)
        : ScalarAssembler<Basis, FEContainer>(basis, fes, dirichFlags) {}

    /** Calculates the vectorial quantity which is requested by fErequirements and returns a reference
     * A zero is written on fixed dofs */
    Eigen::VectorXd &getVector(const RequirementType &fErequirements) { return getVectorImpl(fErequirements); }

    /** Calculates the vectorial quantity which is requested by fErequirements and returns a reference
     * This vector has a reduced size by the number of fixed degrees of freedom */
    Eigen::VectorXd &getReducedVector(const RequirementType &fErequirements) {
      return getReducedVectorImpl(fErequirements);
    }

  private:
    Eigen::VectorXd &getVectorImpl(const RequirementType &fErequirements) {
      vec.setZero(this->size());
      Eigen::VectorXd vecLocal;
      std::vector<GlobalIndex> dofs;
      for (auto &fe : this->finiteElements()) {
        vecLocal.setZero(fe.size());
        dofs.resize(0);
        fe.calculateVector(fErequirements, vecLocal);
        fe.globalIndices(dofs);
        for (int i = 0; auto id : dofs) {
          vec(id[0]) += vecLocal(i);
          ++i;
        }
      }
      for (auto i = 0U; i < this->size(); ++i) {
        if (this->isConstrained(i)) vec[i] = 0;
      }

      return vec;
    }

    Eigen::VectorXd &getReducedVectorImpl(const RequirementType &fErequirements) {
      vecRed.setZero(this->reducedSize());
      int reducedCounter = 0;
      Eigen::VectorXd vecLocal;
      std::vector<GlobalIndex> dofs;
      for (auto &fe : this->finiteElements()) {
        vecLocal.setZero(fe.size());
        dofs.resize(0);
        fe.calculateVector(fErequirements, vecLocal);
        fe.globalIndices(dofs);
        assert(static_cast<long int>(dofs.size()) == vecLocal.size() && "The returned vector has wrong rowSize!");
        for (int i = 0; auto &&dofIndex : dofs) {
          if (this->isConstrained(dofIndex[0])) {
            ++reducedCounter;
            ++i;
            continue;
          } else
            vecRed(dofIndex[0] - this->constraintsBelow(dofIndex[0])) += vecLocal[i++];
        }
      }
      return vecRed;
    }

    Eigen::VectorXd vec{};
    Eigen::VectorXd vecRed{};
  };

  /** SparseFlatAssembler assembles matrix quantities using a flat basis Indexing strategy
   * The matrix is stored in a sparse matrix format. This format is exploited during the assembly process
   */
  template <typename Basis, typename FEContainer>
  class SparseFlatAssembler : public VectorFlatAssembler<Basis, FEContainer> {
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    using GlobalIndex     = typename FEContainer::value_type::GlobalIndex;

  public:
    SparseFlatAssembler(const Basis &basis, const FEContainer &fes, const std::vector<bool> &dirichFlags)
        : VectorFlatAssembler<Basis, FEContainer>(basis, fes, dirichFlags) {}

    using GridView = typename Basis::GridView;

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * A zero is written on fixed dofs rows and columns and a one is written on the diagonal */
    Eigen::SparseMatrix<double> &getMatrix(const RequirementType &fErequirements) {
      return getMatrixImpl(fErequirements);
    }

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * The size of the matrix has the size of the free degrees of freedom */
    Eigen::SparseMatrix<double> &getReducedMatrix(const RequirementType &fErequirements) {
      return getReducedMatrixImpl(fErequirements);
    }

  private:
    Eigen::SparseMatrix<double> &getMatrixImpl(const RequirementType &fErequirements) {
      if (!isOccupationPatternCreated) createOccupationPattern();
      if (!arelinearDofsPerElementCreated) createlinearDofsPerElement();
      spMat.coeffs().setZero();
      Eigen::MatrixXd A;
      for (size_t elementIndex = 0; const auto &fe : this->finiteElements()) {
        A.setZero(fe.size(), fe.size());
        fe.calculateMatrix(fErequirements, A);
        assert(std::sqrt(elementLinearIndices[elementIndex].size()) == A.rows()
               && "The returned matrix has wrong rowSize!");
        assert(std::sqrt(elementLinearIndices[elementIndex].size()) == A.cols()
               && "The returned matrix has wrong colSize!");
        for (Eigen::Index linearIndex = 0; double matrixEntry : A.reshaped())
          spMat.coeffs()(elementLinearIndices[elementIndex][linearIndex++]) += matrixEntry;
        ++elementIndex;
      }
      for (auto i = 0U; i < this->size(); ++i)
        if (this->isConstrained(i)) spMat.col(i) *= 0;
      for (auto i = 0U; i < this->size(); ++i)
        if (this->isConstrained(i)) spMat.row(i) *= 0;
      for (auto i = 0U; i < this->size(); ++i)
        if (this->isConstrained(i)) spMat.diagonal()[i] = 1;
      return spMat;
    }

    Eigen::SparseMatrix<double> &getReducedMatrixImpl(const RequirementType &fErequirements) {
      if (!isReducedOccupationPatternCreated) createReducedOccupationPattern();
      if (!arelinearReducedDofsPerElementCreated) createlinearDofsPerElementReduced();
      spMatReduced.coeffs().setZero();
      Eigen::MatrixXd A;
      std::vector<GlobalIndex> dofs;
      for (size_t elementIndex = 0; const auto &fe : this->finiteElements()) {
        A.setZero(fe.size(), fe.size());
        dofs.resize(0);
        fe.calculateMatrix(fErequirements, A);
        fe.globalIndices(dofs);
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

    /** Calculates the non-zero entries in the full sparse matrix and passed them to the underlying eigen sparse matrix
     *
     * https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
     */
    void createOccupationPattern() {
      spMat.resize(this->size(), this->size());
      std::vector<Eigen::Triplet<double>> vectorOfTriples;

      vectorOfTriples.reserve(this->estimateOfConnectivity());
      std::vector<GlobalIndex> dofs;
      for (auto &&fe : this->finiteElements()) {
        dofs.resize(0);
        fe.globalIndices(dofs);
        for (auto idi : dofs)
          for (auto idj : dofs)
            vectorOfTriples.emplace_back(idi[0], idj[0], 0.0);
      }

      spMat.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
      isOccupationPatternCreated = true;
    }

    /** Calculates the non-zero entries in the sparse matrix and passed them to the underlying eigen sparse matrix
     * The size of the matrix has the size of the free degrees of freedom
     * https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
     */
    void createReducedOccupationPattern() {
      spMatReduced.resize(this->reducedSize(), this->reducedSize());
      std::vector<Eigen::Triplet<double>> vectorOfTriples;
      const int estimateOfConnectivity = 8;
      using std::size;

      vectorOfTriples.reserve(this->estimateOfConnectivity());
      std::vector<GlobalIndex> dofs;
      for (auto &fe : this->finiteElements()) {
        dofs.resize(0);
        fe.globalIndices(dofs);
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

    /** This function save the dof indices of each element in the vector elementLinearIndices */
    void createlinearDofsPerElement() {
      std::vector<GlobalIndex> dofs;
      for (auto &&fe : this->finiteElements()) {
        dofs.resize(0);
        fe.globalIndices(dofs);
        elementLinearIndices.emplace_back(Dune::Power<2>::eval(dofs.size()));
        for (Eigen::Index linearIndexOfElement = 0; auto &&c : dofs)
          for (auto &&r : dofs)
            elementLinearIndices.back()[linearIndexOfElement++] = spMat.getLinearIndex(r[0], c[0]);
      }
      arelinearDofsPerElementCreated = true;
    }

    /** This function save the dof indices of each element in the vector elementLinearIndices but excludes fixed dofs */
    void createlinearDofsPerElementReduced() {
      std::vector<GlobalIndex> dofs;
      for (auto &&fe : this->finiteElements()) {
        dofs.resize(0);
        fe.globalIndices(dofs);
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
      arelinearReducedDofsPerElementCreated = true;
    }

  private:
    Eigen::SparseMatrix<double> spMat;
    Eigen::SparseMatrix<double> spMatReduced;
    bool isOccupationPatternCreated{false};
    bool isReducedOccupationPatternCreated{false};
    bool arelinearDofsPerElementCreated{false};
    bool arelinearReducedDofsPerElementCreated{false};
    std::vector<std::vector<Eigen::Index>> elementLinearIndices;
    std::vector<std::vector<Eigen::Index>> elementLinearReducedIndices;
  };


  /** DenseFlatAssembler assembles matrix quantities using a flat basis Indexing strategy
   * The matrix is stored in a dense matrix format. This format is exploited during the assembly process
   */
  template <typename Basis, typename FEContainer>  // requires Ikarus::Concepts::FlatIndexBasis<BasisEmbedded>
  class DenseFlatAssembler : public VectorFlatAssembler<Basis, FEContainer> {
  public:
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    using GlobalIndex     = typename FEContainer::value_type::GlobalIndex;
    explicit DenseFlatAssembler(const Basis &basis, const FEContainer &fes, const std::vector<bool> &dirichFlags)
        : VectorFlatAssembler<Basis, FEContainer>(basis, fes, dirichFlags) {}



    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * A zero is written on fixed dofs rows and columns and a one is written on the diagonal */
    Eigen::MatrixXd &getMatrix(const RequirementType &fErequirements) { return getMatrixImpl(fErequirements); }

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * The size of the matrix has the size of the free degrees of freedom */
    Eigen::MatrixXd &getReducedMatrix(const RequirementType &fErequirements) {
      return getReducedMatrixImpl(fErequirements);
    }

  private:
    Eigen::MatrixXd &getReducedMatrixImpl(const RequirementType &fErequirements) {
      matRed.setZero(this->reducedSize(), this->reducedSize());
      Eigen::MatrixXd matLocal;
      std::vector<GlobalIndex> dofs;
      for (auto &fe : this->finiteElements()) {
        matLocal.setZero(fe.size(), fe.size());
        dofs.resize(0);
        fe.calculateMatrix(fErequirements, matLocal);
        fe.globalIndices(dofs);
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

    Eigen::MatrixXd &getMatrixImpl(const RequirementType &fErequirements) {
      mat.setZero(this->size(), this->size());
      Eigen::MatrixXd matLocal;
      std::vector<GlobalIndex> dofs;
      for (auto &fe : this->finiteElements()) {
        matLocal.setZero(fe.size(), fe.size());
        dofs.resize(0);
        fe.calculateMatrix(fErequirements, matLocal);
        fe.globalIndices(dofs);
        for (auto i = 0; auto idi : dofs) {
          for (auto j = 0; auto idj : dofs) {
            mat(idi[0], idj[0]) += matLocal(i, j);
            ++j;
          }
          ++i;
        }
      }
      for (auto i = 0U; i < this->size(); ++i)
        if (this->isConstrained(i)) mat.col(i).setZero();
      for (auto i = 0U; i < this->size(); ++i)
        if (this->isConstrained(i)) mat.row(i).setZero();
      for (auto i = 0U; i < this->size(); ++i)
        if (this->isConstrained(i)) mat(i, i) = 1;
      return mat;
    }

    Eigen::MatrixXd mat{};
    Eigen::MatrixXd matRed{};
  };

}  // namespace Ikarus