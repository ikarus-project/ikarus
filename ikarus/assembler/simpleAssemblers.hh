// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <ranges>
#include <utility>

#include <dune/common/math.hh>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus {

  /*!
   * The FlatAssemblerBase takes care of common subtask done by flat assemblers
   */
  template <typename Basis, typename FEContainer>
  class FlatAssemblerBase {
  public:
    using GridView = typename Basis::GridView;
    FlatAssemblerBase(const FEContainer &fes, const DirichletValues<Basis> &p_dirichletValues)
        : feContainer{fes}, dirichletValues{&p_dirichletValues} {
      constraintsBelow_.reserve(dirichletValues->size());
      size_t counter = 0;
      for (auto iv : std::ranges::iota_view{size_t(0), dirichletValues->size()}) {
        constraintsBelow_.emplace_back(counter);
        if (dirichletValues->isConstrained(iv)) ++counter;
      }
      fixedDofs = dirichletValues->fixedDOFsize();
    }

    /**  Returns the size of the free degrees of freeedom, which are not fixed by a dirichlet boundary condition */
    size_t reducedSize() { return dirichletValues->size() - fixedDofs; }

    /**  Returns the size of nodes, i.e. the number of degrees of freedom */
    size_t size() { return dirichletValues->size(); }

    /**  Creates a the fullsized vector of size #Dof and inserts the values of reduced Vector at the "free" degrees
     * of freedom and
     * writes a zero for the fixed doffs */
    Eigen::VectorXd createFullVector(const Eigen::VectorXd &reducedVector);

    /**  Returns the container of finite elements */
    auto &finiteElements() const { return feContainer; }

    /**  Returns the number of constraints below a given degrees of freedom index */
    size_t constraintsBelow(size_t i) const { return constraintsBelow_[i]; }

    /**  Returns the boolean if a given degree of freedom is fixed */
    bool isConstrained(size_t i) const { return dirichletValues->isConstrained(i); }

    /**  Coarse estimate of node connectivity, i.e. this relates the bandwidth of an sparse matrix.
     * This estimate is designed that it overestimates the real connectivity since it should
     * only be used for allocating vectors */
    size_t estimateOfConnectivity() const { return dirichletValues->basis().gridView().size(GridView::dimension) * 8; }

  private:
    FEContainer const &feContainer;
    DirichletValues<Basis> const *dirichletValues;
    std::vector<size_t> constraintsBelow_{};
    size_t fixedDofs{};
  };

  /** ScalarAssembler assembles scalar quantities */
  template <typename Basis, typename FEContainer>
  class ScalarAssembler : public FlatAssemblerBase<Basis, FEContainer> {
    using RequirementType = typename FEContainer::value_type::FERequirementType;

  public:
    ScalarAssembler(const FEContainer &fes, const DirichletValues<Basis> &dirichletValues_)
        : FlatAssemblerBase<Basis, FEContainer>(fes, dirichletValues_) {}

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
    VectorFlatAssembler(const FEContainer &fes, const DirichletValues<Basis> &dirichletValues_)
        : ScalarAssembler<Basis, FEContainer>(fes, dirichletValues_) {}

    /** Calculates the vectorial quantity which is requested by fErequirements and returns a reference
     * A zero is written on fixed dofs */
    Eigen::VectorXd &getVector(const RequirementType &fErequirements) { return getVectorImpl(fErequirements); }

    /** Calculates the vectorial quantity which is requested by fErequirements and returns a reference
     * This vector has a reduced size by the number of fixed degrees of freedom */
    Eigen::VectorXd &getReducedVector(const RequirementType &fErequirements) {
      return getReducedVectorImpl(fErequirements);
    }

  private:
    Eigen::VectorXd &getVectorImpl(const RequirementType &fErequirements);
    Eigen::VectorXd &getReducedVectorImpl(const RequirementType &fErequirements);

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
    SparseFlatAssembler(const FEContainer &fes, const DirichletValues<Basis> &dirichletValues_)
        : VectorFlatAssembler<Basis, FEContainer>(fes, dirichletValues_) {}

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
    Eigen::SparseMatrix<double> &getMatrixImpl(const RequirementType &fErequirements);
    Eigen::SparseMatrix<double> &getReducedMatrixImpl(const RequirementType &fErequirements);

    /** Calculates the non-zero entries in the full sparse matrix and passed them to the underlying eigen sparse matrix
     * https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
     */
    void createOccupationPattern();

    /** Calculates the non-zero entries in the sparse matrix and passed them to the underlying eigen sparse matrix
     * The size of the matrix has the size of the free degrees of freedom
     * https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
     */
    void createReducedOccupationPattern();

    /** This function save the dof indices of each element in the vector elementLinearIndices */
    void createlinearDofsPerElement();

    /** This function save the dof indices of each element in the vector elementLinearIndices but excludes fixed dofs */
    void createlinearDofsPerElementReduced();

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
    explicit DenseFlatAssembler(const FEContainer &fes, const DirichletValues<Basis> &dirichletValues_)
        : VectorFlatAssembler<Basis, FEContainer>(fes, dirichletValues_) {}

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * A zero is written on fixed dofs rows and columns and a one is written on the diagonal */
    Eigen::MatrixXd &getMatrix(const RequirementType &fErequirements) { return getMatrixImpl(fErequirements); }

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * The size of the matrix has the size of the free degrees of freedom */
    Eigen::MatrixXd &getReducedMatrix(const RequirementType &fErequirements) {
      return getReducedMatrixImpl(fErequirements);
    }

  private:
    Eigen::MatrixXd &getReducedMatrixImpl(const RequirementType &fErequirements);
    Eigen::MatrixXd &getMatrixImpl(const RequirementType &fErequirements);

    Eigen::MatrixXd mat{};
    Eigen::MatrixXd matRed{};
  };

}  // namespace Ikarus

#include "simpleAssemblers.inl"
