// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <ranges>
#include <utility>

#include <dune/common/math.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus {

  /*!
   * The FlatAssemblerBase takes care of common subtask done by flat assemblers
   */
  template <typename FEContainer_, typename DirichletValuesType_>
  class FlatAssemblerBase {
  public:
    using FEContainerRaw    = std::remove_cvref_t<FEContainer_>;
    using FERequirementType = typename FEContainerRaw::value_type::FERequirementType;
    using GlobalIndex       = typename FEContainerRaw::value_type::GlobalIndex;
    using Basis             = typename DirichletValuesType_::Basis;
    using GridView          = typename Basis::GridView;
    using FEContainer       = FEContainer_;
    // If we get a reference we store the reference but if we get an r-value we store by value
    using FEContainerType     = std::conditional_t<std::is_reference_v<FEContainer>, const FEContainer, FEContainer>;
    using DirichletValuesType = DirichletValuesType_;
    FlatAssemblerBase(FEContainer &&fes, const DirichletValuesType &p_dirichletValues)
        : feContainer{std::forward<FEContainer>(fes)}, dirichletValues{&p_dirichletValues} {
      constraintsBelow_.reserve(dirichletValues->size());
      size_t counter = 0;
      for (auto iv : std::ranges::iota_view{decltype(dirichletValues->size())(0), dirichletValues->size()}) {
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
    Eigen::VectorXd createFullVector(Eigen::Ref<const Eigen::VectorXd> reducedVector);

    /**  Returns the container of finite elements */
    auto &finiteElements() const { return feContainer; }

    /**  Returns the number of constraints below a given degrees of freedom index */
    [[nodiscard]] size_t constraintsBelow(size_t i) const { return constraintsBelow_[i]; }

    /**  Returns the boolean if a given degree of freedom is fixed */
    [[nodiscard]] bool isConstrained(size_t i) const { return dirichletValues->isConstrained(i); }

    /**  Coarse estimate of node connectivity, i.e. this relates the bandwidth of an sparse matrix.
     * This estimate is designed that it overestimates the real connectivity since it should
     * only be used for allocating vectors */
    [[nodiscard]] size_t estimateOfConnectivity() const {
      return dirichletValues->basis().gridView().size(GridView::dimension) * 8;
    }

  private:
    FEContainerType feContainer;
    DirichletValuesType const *dirichletValues;
    std::vector<size_t> constraintsBelow_{};
    size_t fixedDofs{};
  };

  template <class T, class DirichletValuesType>
  FlatAssemblerBase(T &&fes, const DirichletValuesType &dirichletValues_) -> FlatAssemblerBase<T, DirichletValuesType>;

  /** ScalarAssembler assembles scalar quantities */
  template <typename FEContainer_, typename DirichletValuesType_>
  class ScalarAssembler : public FlatAssemblerBase<FEContainer_, DirichletValuesType_> {
    using FEContainerRaw = std::remove_cvref_t<FEContainer_>;

    using FERequirementType   = typename FEContainerRaw::value_type::FERequirementType;
    using GlobalIndex         = typename FEContainerRaw::value_type::GlobalIndex;
    using Basis               = typename DirichletValuesType_::Basis;
    using FEContainer         = FEContainer_;
    using DirichletValuesType = DirichletValuesType_;

  public:
    ScalarAssembler(FEContainer &&fes, const DirichletValuesType &dirichletValues_)
        : FlatAssemblerBase<FEContainer, DirichletValuesType>(std::forward<FEContainer>(fes), dirichletValues_) {}

    /** Calculates the scalar quantity which is requested by fErequirements and returns a reference */
    double &getScalar(const FERequirementType &fErequirements) { return getScalarImpl(fErequirements); }

  private:
    double &getScalarImpl(const FERequirementType &fErequirements) {
      scal = 0.0;
      for (auto &fe : this->finiteElements()) {
        scal += fe.calculateScalar(fErequirements);
      }
      return scal;
    }

    double scal{0.0};
  };

  template <class T, class DirichletValuesType>
  ScalarAssembler(T &&fes, const DirichletValuesType &dirichletValues_) -> ScalarAssembler<T, DirichletValuesType>;

  /** VectorFlatAssembler assembles vector quantities using a flat basis Indexing strategy */
  template <typename FEContainer_, typename DirichletValuesType_>
  class VectorFlatAssembler : public ScalarAssembler<FEContainer_, DirichletValuesType_> {
    using FEContainerRaw = std::remove_cvref_t<FEContainer_>;

    using FERequirementType   = typename FEContainerRaw::value_type::FERequirementType;
    using GlobalIndex         = typename FEContainerRaw::value_type::GlobalIndex;
    using Basis               = typename DirichletValuesType_::Basis;
    using FEContainer         = FEContainer_;
    using DirichletValuesType = DirichletValuesType_;

  public:
    VectorFlatAssembler(FEContainer &&fes, const DirichletValuesType &dirichletValues_)
        : ScalarAssembler<FEContainer, DirichletValuesType>(std::forward<FEContainer>(fes), dirichletValues_) {}

    /** Calculates the vectorial quantity which is requested by fErequirements and returns a reference
     * A zero is written on fixed dofs */
    Eigen::VectorXd &getVector(const FERequirementType &fErequirements) { return getVectorImpl(fErequirements); }

    /** Calculates the vectorial quantity which is requested by fErequirements and returns a reference
     * This vector has a reduced size by the number of fixed degrees of freedom */
    Eigen::VectorXd &getReducedVector(const FERequirementType &fErequirements) {
      return getReducedVectorImpl(fErequirements);
    }

  private:
    Eigen::VectorXd &getVectorImpl(const FERequirementType &fErequirements);
    Eigen::VectorXd &getReducedVectorImpl(const FERequirementType &fErequirements);

    Eigen::VectorXd vec{};
    Eigen::VectorXd vecRed{};
  };

  template <class T, class DirichletValuesType>
  VectorFlatAssembler(T &&fes, const DirichletValuesType &dirichletValues_)
      -> VectorFlatAssembler<T, DirichletValuesType>;

  /** SparseFlatAssembler assembles matrix quantities using a flat basis Indexing strategy
   * The matrix is stored in a sparse matrix format. This format is exploited during the assembly process
   */
  template <typename FEContainer_, typename DirichletValuesType_>
  class SparseFlatAssembler : public VectorFlatAssembler<FEContainer_, DirichletValuesType_> {
  public:
    using FEContainerRaw = std::remove_cvref_t<FEContainer_>;

    using FERequirementType   = typename FEContainerRaw::value_type::FERequirementType;
    using GlobalIndex         = typename FEContainerRaw::value_type::GlobalIndex;
    using Basis               = typename DirichletValuesType_::Basis;
    using FEContainer         = FEContainer_;
    using DirichletValuesType = DirichletValuesType_;

    SparseFlatAssembler(FEContainer &&fes, const DirichletValuesType &dirichletValues_)
        : VectorFlatAssembler<FEContainer, DirichletValuesType>(std::forward<FEContainer>(fes), dirichletValues_) {}

    using GridView = typename Basis::GridView;

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * A zero is written on fixed dofs rows and columns and a one is written on the diagonal */
    Eigen::SparseMatrix<double> &getMatrix(const FERequirementType &fErequirements) {
      return getMatrixImpl(fErequirements);
    }

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * The size of the matrix has the size of the free degrees of freedom */
    Eigen::SparseMatrix<double> &getReducedMatrix(const FERequirementType &fErequirements) {
      return getReducedMatrixImpl(fErequirements);
    }

  private:
    Eigen::SparseMatrix<double> &getMatrixImpl(const FERequirementType &fErequirements);
    Eigen::SparseMatrix<double> &getReducedMatrixImpl(const FERequirementType &fErequirements);

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

  template <class T, class DirichletValuesType>
  SparseFlatAssembler(T &&fes, const DirichletValuesType &dirichletValues_)
      -> SparseFlatAssembler<T, DirichletValuesType>;

  /** DenseFlatAssembler assembles matrix quantities using a flat basis Indexing strategy
   * The matrix is stored in a dense matrix format. This format is exploited during the assembly process
   */
  template <typename FEContainer_,
            typename DirichletValuesType_>  // requires Ikarus::Concepts::FlatIndexBasis<BasisEmbedded>
  class DenseFlatAssembler : public VectorFlatAssembler<FEContainer_, DirichletValuesType_> {
  public:
    using FEContainerRaw = std::remove_cvref_t<FEContainer_>;

    using FERequirementType   = typename FEContainerRaw::value_type::FERequirementType;
    using GlobalIndex         = typename FEContainerRaw::value_type::GlobalIndex;
    using Basis               = typename DirichletValuesType_::Basis;
    using FEContainer         = FEContainer_;
    using DirichletValuesType = DirichletValuesType_;

    explicit DenseFlatAssembler(FEContainer &&fes, const DirichletValuesType &dirichletValues_)
        : VectorFlatAssembler<FEContainer, DirichletValuesType>(std::forward<FEContainer>(fes), dirichletValues_) {}

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * A zero is written on fixed dofs rows and columns and a one is written on the diagonal */
    Eigen::MatrixXd &getMatrix(const FERequirementType &fErequirements) { return getMatrixImpl(fErequirements); }

    /** Calculates the matrix quantity which is requested by fErequirements and returns a reference
     * The size of the matrix has the size of the free degrees of freedom */
    Eigen::MatrixXd &getReducedMatrix(const FERequirementType &fErequirements) {
      return getReducedMatrixImpl(fErequirements);
    }

  private:
    Eigen::MatrixXd &getReducedMatrixImpl(const FERequirementType &fErequirements);
    Eigen::MatrixXd &getMatrixImpl(const FERequirementType &fErequirements);

    Eigen::MatrixXd mat{};
    Eigen::MatrixXd matRed{};
  };

  // https://en.cppreference.com/w/cpp/language/class_template_argument_deduction
  template <class T, class DirichletValuesType>
  DenseFlatAssembler(T &&fes, const DirichletValuesType &dirichletValues_)
      -> DenseFlatAssembler<T, DirichletValuesType>;

}  // namespace Ikarus

#include "simpleAssemblers.inl"
