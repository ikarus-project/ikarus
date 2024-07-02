// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file simpleassemblers.hh
 * \brief Defines several assemblers for finite element assembly
 */

#pragma once

#include <ranges>
#include <utility>

#include <dune/functions/backends/istlvectorbackend.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/assembler/interface.hh>
#include <ikarus/finiteelements/fehelper.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/dirichletvalues.hh>

namespace Ikarus {
/**
 * \class FlatAssemblerBase
 * \brief The FlatAssemblerBase takes care of common subtasks done by flat assemblers.
 * \ingroup assembler
 *
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class FlatAssemblerBase
{
public:
  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using FERequirement  = typename FEContainerRaw::value_type::Requirement;
  ///< Type of the finite element requirement.
  using GlobalIndex     = typename FEContainerRaw::value_type::GlobalIndex; ///< Type of the global index.
  using Basis           = typename DV::Basis;                               ///< Type of the basis.
  using GridView        = typename Basis::GridView;                         ///< Type of the grid view.
  using FEContainer     = FEC;                                              ///< Type of the finite element container.
  using FEContainerType = std::conditional_t<std::is_reference_v<FEContainer>, const FEContainer, FEContainer>;
  ///< Type of the finite element container (reference or by value).
  using DirichletValuesType = DV; ///< Type of the Dirichlet values.
  using SizeType =
      typename DirichletValuesType::FlagsType::size_type; ///< size_type of the container storing Dirichlet flags

  /**
   * \brief Constructor for FlatAssemblerBase.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  FlatAssemblerBase(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : feContainer_{std::forward<FEContainer>(fes)},
        dirichletValues_{&dirichletValues} {
    constraintsBelow_.reserve(dirichletValues_->size());
    size_t counter = 0;
    for (auto iv : std::ranges::iota_view{decltype(dirichletValues_->size())(0), dirichletValues_->size()}) {
      constraintsBelow_.emplace_back(counter);
      if (dirichletValues_->isConstrained(iv))
        ++counter;
    }
    fixedDofs_ = dirichletValues_->fixedDOFsize();
  }

  /**
   * \brief Returns the size of the free degrees of freedom, which are not fixed by a Dirichlet boundary condition.
   * \return Size of the reduced degrees of freedom.
   */
  size_t reducedSize() { return dirichletValues_->size() - fixedDofs_; }

  /**
   * \brief Returns the size of nodes, i.e., the number of degrees of freedom.
   * \return Size of the degrees of freedom.
   */
  size_t size() { return dirichletValues_->size(); }

  /**
   * \brief Creates the full-sized vector of size Dof and inserts the values of a reduced vector at the "free"
   * degrees of freedom and writes a zero for the fixed doffs.
   *
   * \param reducedVector Reference to the reduced vector.
   * \return Eigen::VectorXd The full-sized vector.
   */
  Eigen::VectorXd createFullVector(Eigen::Ref<const Eigen::VectorXd> reducedVector);

  /**
   * \brief Returns the container of finite elements.
   * \return Reference to the finite element container.
   */
  auto& finiteElements() const { return feContainer_; }

  /**
   * \brief Returns the dirichlet value object.
   * \return Reference to the dirichlet value object.
   */
  const auto& dirichletValues() const { return *dirichletValues_; }

  /**
   * \brief Returns the number of constraints below a given degrees of freedom index.
   *
   * \param i Index of the degree of freedom.
   * \return Number of constraints below the given index.
   */
  [[nodiscard]] size_t constraintsBelow(SizeType i) const { return constraintsBelow_[i]; }

  /**
   * \brief Returns true if a given degree of freedom is fixed by a Dirichlet boundary condition.
   *
   * \param i Index of the degree of freedom.
   * \return True if the degree of freedom is fixed; false otherwise.
   */
  [[nodiscard]] bool isConstrained(SizeType i) const { return dirichletValues_->isConstrained(i); }

  /**
   * \brief Coarse estimate of node connectivity, i.e., this relates to the bandwidth of a sparse matrix.
   * This estimate overestimates the real connectivity and should only be used for allocating vectors.
   *
   * \return Size_t Coarse estimate of node connectivity.
   */
  [[nodiscard]] size_t estimateOfConnectivity() const {
    return dirichletValues_->basis().gridView().size(GridView::dimension) * 8;
  }

  using AffordanceCollectionType = AffordanceCollection<ScalarAffordance, VectorAffordance, MatrixAffordance>;

  /**
   * \brief Binds the assembler to a set of finite element requirement and affordance.
   *
   * \param req Reference to the finite element requirement.
   * \param affordanceCollection The affordance collection
   */
  void bind(const FERequirement& req, AffordanceCollectionType affordanceCollection,
            DBCOption dbcOption = DBCOption::Full) {
    req_         = std::make_optional<FERequirement>(req);
    affordances_ = std::make_optional<AffordanceCollectionType>(affordanceCollection);
    dBCOption_   = std::make_optional<DBCOption>(dbcOption);
  }

  /**
   * \brief Binds the assembler to a  finite element requirement.
   *
   * \param req Reference to the finite element requirement.
   */
  void bind(const FERequirement& req) { req_ = std::make_optional<FERequirement>(req); }

  /**
   * \brief Binds the assembler to an affordance collection.
   *
   * \param affordanceCollection The affordance collection
   */
  void bind(AffordanceCollectionType affordanceCollection) {
    affordances_ = std::make_optional<AffordanceCollectionType>(affordanceCollection);
  }

  /**
   * \brief Binds the assembler to an affordance collection.
   *
   * \param dbcOption The EnforcingDBC option
   */
  void bind(DBCOption dbcOption) { dBCOption_ = std::make_optional<DBCOption>(dbcOption); }

  /**
   * \brief Returns true if the assembler is bound to a finite element requirement and affordance.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool bound() const {
    if (boundToRequirement() and boundToAffordanceCollection() and boundToDBCOption())
      return true;
    else
      DUNE_THROW(Dune::InvalidStateException, "The assembler is not bound to a requirement, affordance or dBCOption.");
  }

  /**
   * \brief Returns true if the assembler is bound to a finite element requirement.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool boundToRequirement() const {
    return req_.has_value();
  }

  /**
   * \brief Returns true if the assembler is bound to an affordance collection.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool boundToAffordanceCollection() const {
    return affordances_.has_value();
  }

  /**
   * \brief Returns true if the assembler is bound to an affordance collection.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool boundToDBCOption() const {
    return dBCOption_.has_value();
  }

  /**
   * \brief Returns the requirement.
   *
   */
  FERequirement& requirement() {
    if (req_.has_value())
      return req_.value();
    else
      DUNE_THROW(Dune::InvalidStateException, "The requirement can only be obtained after binding");
  }

  /**
   * \brief Returns the affordance.
   *
   */
  AffordanceCollectionType affordanceCollection() const {
    if (affordances_.has_value())
      return affordances_.value();
    else
      DUNE_THROW(Dune::InvalidStateException, "The affordance can only be obtained after binding");
  }

  /**
   * \brief Returns the dirichlet boundary condition enforcement option.
   *
   */
  DBCOption dBCOption() const {
    if (dBCOption_.has_value())
      return dBCOption_.value();
    else
      DUNE_THROW(Dune::InvalidStateException, "The dBCOption can only be obtained after binding");
  }

private:
  FEContainerType feContainer_;
  const DirichletValuesType* dirichletValues_;
  std::optional<FERequirement> req_;
  std::optional<AffordanceCollectionType> affordances_;
  std::vector<size_t> constraintsBelow_{};
  size_t fixedDofs_{};
  std::optional<DBCOption> dBCOption_;
};

#ifndef DOXYGEN
template <class T, class DirichletValuesType>
FlatAssemblerBase(T&& fes, const DirichletValuesType& dirichletValues) -> FlatAssemblerBase<T, DirichletValuesType>;
#endif

/**
 * \class ScalarFlatAssembler
 * \brief ScalarFlatAssembler assembles scalar quantities.
 * \ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class ScalarFlatAssembler : public ScalarAssembler<ScalarFlatAssembler<FEC, DV>, FEC, DV, double>,
                            public FlatAssemblerBase<FEC, DV>
{
protected:
  using Base = FlatAssemblerBase<FEC, DV>; ///< Type alias for the base class.
  friend ScalarAssembler<ScalarFlatAssembler, FEC, DV, double>;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;
  using typename ScalarAssembler<ScalarFlatAssembler, FEC, DV, double>::ScalarType;

  /**
   * \brief Constructor for ScalarFlatAssembler.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  ScalarFlatAssembler(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : FlatAssemblerBase<FEC, DV>(std::forward<FEContainer>(fes), dirichletValues) {}

protected:
  ScalarType& getScalarImpl(const FERequirement& feRequirements, ScalarAffordance affordance);

  ScalarType scal_{0.0};
};

#ifndef DOXYGEN
template <class T, class DirichletValuesType>
ScalarFlatAssembler(T&& fes, const DirichletValuesType& dirichletValues) -> ScalarFlatAssembler<T, DirichletValuesType>;
#endif

/**
 * \class VectorFlatAssembler
 * \brief VectorFlatAssembler assembles vector quantities using a flat basis Indexing strategy.
 *\ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class VectorFlatAssembler : public ScalarFlatAssembler<FEC, DV>,
                            public VectorAssembler<VectorFlatAssembler<FEC, DV>, FEC, DV, Eigen::VectorXd>
{
protected:
  using Base = ScalarFlatAssembler<FEC, DV>; ///< Type alias for the base class.
  friend VectorAssembler<VectorFlatAssembler, FEC, DV, Eigen::VectorXd>;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;

  using typename Base::ScalarType;
  using typename VectorAssembler<VectorFlatAssembler, FEC, DV, Eigen::VectorXd>::VectorType;

  /**
   * \brief Constructor for VectorFlatAssembler.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  VectorFlatAssembler(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : Base(std::forward<FEContainer>(fes), dirichletValues) {}

protected:
  void assembleRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance, VectorType& assemblyVec);
  VectorType& getRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);
  VectorType& getVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);

  VectorType& getReducedVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);

  VectorType vecRaw_{}; ///< Raw vector without changes for dirichlet degrees of freedom
  VectorType vec_{};    ///< Vector quantity.
  VectorType vecRed_{}; ///< Reduced vector quantity.
};

#ifndef DOXYGEN
template <class T, class DirichletValuesType>
VectorFlatAssembler(T&& fes, const DirichletValuesType& dirichletValues) -> VectorFlatAssembler<T, DirichletValuesType>;
#endif

/**
 * \class SparseFlatAssembler
 * \brief SparseFlatAssembler assembles matrix quantities using a flat basis Indexing strategy.
 * The matrix is stored in a sparse matrix format. This format is exploited during the assembly process.
 * \ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class SparseFlatAssembler : public MatrixAssembler<SparseFlatAssembler<FEC, DV>, FEC, DV, Eigen::SparseMatrix<double>>,
                            public VectorFlatAssembler<FEC, DV>
{
protected:
  using Base = VectorFlatAssembler<FEC, DV>; ///< Type alias for the base class.
  friend MatrixAssembler<SparseFlatAssembler, FEC, DV, Eigen::SparseMatrix<double>>;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;
  using typename Base::ScalarType;
  using typename Base::VectorType;
  using typename MatrixAssembler<SparseFlatAssembler, FEC, DV, Eigen::SparseMatrix<double>>::MatrixType;

  /**
   * \brief Constructor for SparseFlatAssembler.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  SparseFlatAssembler(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : Base(std::forward<FEContainer>(fes), dirichletValues) {}

private:
  void assembleRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance, MatrixType& assemblyMat);

protected:
  MatrixType& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  MatrixType& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  MatrixType& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);

private:
  /** Calculates the non-zero entries in the full sparse matrix and passes them to the underlying Eigen sparse matrix.
   */
  void createOccupationPattern(MatrixType& assemblyMat);

  /** Calculates the non-zero entries in the sparse matrix and passes them to the underlying Eigen sparse matrix.
   * The size of the matrix has the size of the free degrees of freedom. */
  void createReducedOccupationPattern(MatrixType& assemblyMat);

  /** Saves the degree of freedom indices of each element in the vector elementLinearIndices. */
  void createLinearDOFsPerElement(MatrixType& assemblyMat);

  /** Saves the degree of freedom indices of each element in the vector elementLinearIndices but excludes fixed
   * degrees of freedom. */
  void createLinearDOFsPerElementReduced(MatrixType& assemblyMat);

  /** Pre-processes the raw sparse matrix before assembly. */
  void preProcessSparseMatrix(MatrixType& assemblyMat);

  /** Pre-processes the reduced sparse matrix before assembly. */
  void preProcessSparseMatrixReduced(MatrixType& assemblyMat);

  MatrixType spMatRaw_;     ///< Raw sparse matrix without changes for dirichlet degrees of freedom.
  MatrixType spMat_;        ///< Sparse matrix.
  MatrixType spMatReduced_; ///< Reduced sparse matrix.
  std::vector<std::vector<Eigen::Index>>
      elementLinearIndices_; ///< Vector storing indices of matrix entries in linear storage
  std::vector<std::vector<Eigen::Index>>
      elementLinearReducedIndices_; ///< Vector storing indices of matrix entries in linear storage
  bool sparsePreProcessorRaw_{}, sparsePreProcessor_{},
      sparsePreProcessorReduced_{}; ///< flags that store if the sparsity pattern construction happened
};

#ifndef DOXYGEN
template <class FEC, class DV>
SparseFlatAssembler(FEC&& fes, const DV& dirichletValues) -> SparseFlatAssembler<FEC, DV>;
#endif

template <typename FEC, typename DV>
auto makeSparseFlatAssembler(FEC&& fes, const DV& dirichletValues) {
  return std::make_shared<SparseFlatAssembler<FEC, DV>>(std::forward<FEC>(fes), dirichletValues);
}

/**
 * \class DenseFlatAssembler
 * \brief DenseFlatAssembler assembles matrix quantities using a flat basis Indexing strategy.
 * The matrix is stored in a dense matrix format. This format is exploited during the assembly process.
 * \ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 * \note Requires Ikarus::Concepts::FlatIndexBasis<BasisEmbedded>.
 */
template <typename FEC, typename DV>
class DenseFlatAssembler : public MatrixAssembler<DenseFlatAssembler<FEC, DV>, FEC, DV, Eigen::MatrixXd>,
                           public VectorFlatAssembler<FEC, DV>
{
protected:
  using Base = VectorFlatAssembler<FEC, DV>; ///< Type alias for the base class.
  friend MatrixAssembler<DenseFlatAssembler, FEC, DV, Eigen::MatrixXd>;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;
  using typename Base::ScalarType;
  using typename Base::VectorType;
  using typename MatrixAssembler<DenseFlatAssembler, FEC, DV, Eigen::MatrixXd>::MatrixType;

  /**
   * \brief Constructor for DenseFlatAssembler.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  explicit DenseFlatAssembler(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : Base(std::forward<FEContainer>(fes), dirichletValues) {}

private:
  void assembleRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance, MatrixType& assemblyMat);

protected:
  MatrixType& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  MatrixType& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  MatrixType& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);

private:
  MatrixType matRaw_{}; ///< Raw dense matrix for assembly.
  MatrixType mat_{};    ///< Dense matrix quantity.
  MatrixType matRed_{}; ///< Reduced dense matrix quantity.
};

#ifndef DOXYGEN
// https://en.cppreference.com/w/cpp/language/class_template_argument_deduction
template <class FEC, class DV>
DenseFlatAssembler(FEC&& fes, const DV& dirichletValues) -> DenseFlatAssembler<FEC, DV>;
#endif

template <typename FEC, typename DV>
auto makeDenseFlatAssembler(FEC&& fes, const DV& dirichletValues) {
  return std::make_shared<DenseFlatAssembler<FEC, DV>>(std::forward<FEC>(fes), dirichletValues);
}
} // namespace Ikarus

#include "simpleassemblers.inl"
