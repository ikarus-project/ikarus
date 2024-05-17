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
   * \brief Returns the number of constraints below a given degrees of freedom index.
   *
   * \param i Index of the degree of freedom.
   * \return Number of constraints below the given index.
   */
  [[nodiscard]] size_t constraintsBelow(size_t i) const { return constraintsBelow_[i]; }

  /**
   * \brief Returns true if a given degree of freedom is fixed by a Dirichlet boundary condition.
   *
   * \param i Index of the degree of freedom.
   * \return True if the degree of freedom is fixed; false otherwise.
   */
  [[nodiscard]] bool isConstrained(size_t i) const { return dirichletValues_->isConstrained(i); }

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
   * \param affo The affordance
   */
  void bind(const FERequirement& req, AffordanceCollectionType affo, EnforcingDBCOption qt = EnforcingDBCOption::Full) {
    req_                = std::make_optional<FERequirement>(req);
    affos_              = std::make_optional<AffordanceCollectionType>(affo);
    enforcingDBCOption_ = std::make_optional<EnforcingDBCOption>(qt);
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
   * \param affo The affordance collection
   */
  void bind(AffordanceCollectionType affo) { affos_ = std::make_optional<AffordanceCollectionType>(affo); }

  /**
   * \brief Binds the assembler to an affordance collection.
   *
   * \param affo The affordance collection
   */
  void bind(EnforcingDBCOption qt) { enforcingDBCOption_ = std::make_optional<EnforcingDBCOption>(qt); }

  /**
   * \brief Returns true if the assembler is bound to a finite element requirement and affordance.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool bound() const {
    return req_.has_value() and affos_.has_value() and enforcingDBCOption_.has_value();
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
    return affos_.has_value();
  }

  /**
   * \brief Returns true if the assembler is bound to an affordance collection.
   *
   * \return True if the assembler is bound; false otherwise.
   */
  [[nodiscard]]
  bool boundToEnforcingDBCOption() const {
    return enforcingDBCOption_.has_value();
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
    if (affos_.has_value())
      return affos_.value();
    else
      DUNE_THROW(Dune::InvalidStateException, "The affordance can only be obtained after binding");
  }

  /**
   * \brief Returns the dirichlet boundary condition enforcement option.
   *
   */
  EnforcingDBCOption enforcingDBCOption() const {
    if (enforcingDBCOption_.has_value())
      return enforcingDBCOption_.value();
    else
      DUNE_THROW(Dune::InvalidStateException, "The enforcingDBCOption can only be obtained after binding");
  }

private:
  FEContainerType feContainer_;
  const DirichletValuesType* dirichletValues_;
  std::optional<FERequirement> req_;
  std::optional<AffordanceCollectionType> affos_;
  std::vector<size_t> constraintsBelow_{};
  size_t fixedDofs_{};
  std::optional<EnforcingDBCOption> enforcingDBCOption_;
};

#ifndef DOXYGEN
template <class T, class DirichletValuesType>
FlatAssemblerBase(T&& fes, const DirichletValuesType& dirichletValues) -> FlatAssemblerBase<T, DirichletValuesType>;
#endif

/**
 * \class ScalarAssembler
 * \brief ScalarAssembler assembles scalar quantities.
 * \ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class ScalarAssembler : public FlatAssemblerBase<FEC, DV>
{
  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using Base           = FlatAssemblerBase<FEC, DV>;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;

  /**
   * \brief Constructor for ScalarAssembler.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  ScalarAssembler(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : FlatAssemblerBase<FEContainer, DirichletValuesType>(std::forward<FEContainer>(fes), dirichletValues) {}

  /**
   * \brief Calculates the scalar quantity requested by feRequirements and affordance.
   *
   * \param feRequirements Reference to the finite element requirements.
   * \param affordance The scalar affordance
   * \return Const reference to the calculated scalar quantity.
   */
  const double& scalar(const FERequirement& feRequirements, ScalarAffordance affordance) {
    return getScalarImpl(feRequirements, affordance);
  }

  /**
   * \brief Calculates the scalar quantity requested by the bound feRequirements and returns a reference.
   *
   * \return Const reference to the calculated scalar quantity.
   */
  const double& scalar() { return getScalarImpl(this->requirement(), this->affordanceCollection().scalarAffordance()); }

private:
  /**
   * \brief Helper function to calculate the scalar quantity based on finite element requirements.
   *
   * \param feRequirements Reference to the finite element requirements.
   * \return Reference to the calculated scalar quantity.
   */
  double& getScalarImpl(const FERequirement& feRequirements, ScalarAffordance affordance) {
    scal_ = 0.0;
    for (auto& fe : this->finiteElements()) {
      scal_ += calculateScalar(fe, feRequirements, affordance);
    }
    return scal_;
  }

  double scal_{0.0};
};

#ifndef DOXYGEN
template <class T, class DirichletValuesType>
ScalarAssembler(T&& fes, const DirichletValuesType& dirichletValues) -> ScalarAssembler<T, DirichletValuesType>;
#endif

/**
 * \class VectorFlatAssembler
 * \brief VectorFlatAssembler assembles vector quantities using a flat basis Indexing strategy.
 *\ingroup assembler
 * \tparam FEC Type of the finite element container.
 * \tparam DV Type of the Dirichlet values.
 */
template <typename FEC, typename DV>
class VectorFlatAssembler : public ScalarAssembler<FEC, DV>
{
  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using Base           = ScalarAssembler<FEC, DV>;

public:
  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;

public:
  /**
   * \brief Constructor for VectorFlatAssembler.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  VectorFlatAssembler(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : ScalarAssembler<FEContainer, DirichletValuesType>(std::forward<FEContainer>(fes), dirichletValues) {}

  /**
   * \brief Calculates the vectorial quantity requested by the  feRequirements and the affordance.
   Depending on the requested EnforcingDBCOption, the raw, reduced or full vector is returned.
    Raw means the degrees of freedom associated with dirichlet boundary conditions are not changed.
    Full means that degrees of freedom associated with dirichlet boundary conditions are set to zero in the vector.
    Reduced means that degrees of freedom associated with dirichlet boundary conditions are removed and the returned
   vector has reduced size.
   *
   * \param feRequirements Reference to the finite element requirements.
   * \param affordance The vector affordance
   * \param qt The EnforcingDBCOption
   * \return Const reference to the calculated vectorial quantity.
   */
  const Eigen::VectorXd& vector(const FERequirement& feRequirements, VectorAffordance affordance,
                                EnforcingDBCOption qt = EnforcingDBCOption::Full) {
    if (qt == EnforcingDBCOption::Raw) {
      return getRawVectorImpl(feRequirements, affordance);
    } else if (qt == EnforcingDBCOption::Reduced) {
      return getReducedVectorImpl(feRequirements, affordance);
    } else if (qt == EnforcingDBCOption::Full) {
      return getVectorImpl(feRequirements, affordance);
    }
    __builtin_unreachable();
  }

  /**
 * \brief Calculates the vectorial quantity requested by the bound feRequirements and the affordance.
 Depending on the requested EnforcingDBCOption, the raw, reduced or full vector is returned.
  Raw means the degrees of freedom associated with dirichlet boundary conditions are not changed.
  Full means that degrees of freedom associated with dirichlet boundary conditions are set to zero in the vector.
  Reduced means that degrees of freedom associated with dirichlet boundary conditions are removed and the returned
 vector has reduced size.
  * \see const Eigen::VectorXd& getVector(const FERequirement& feRequirements, VectorAffordance affordance,
 EnforcingDBCOption qt)
 * \param qt The EnforcingDBCOption
 * \return Const reference to the calculated vectorial quantity.
 */
  const Eigen::VectorXd& vector(EnforcingDBCOption qt) {
    return vector(this->requirement(), this->affordanceCollection().vectorAffordance(), qt);
  }

  /**
* \brief Calculates the vectorial quantity requested by the bound feRequirements,  the affordance and the
enforcingDBCOption. Depending on the EnforcingDBCOption, the raw, reduced or full vector is returned. Raw
means the degrees of freedom associated with dirichlet boundary conditions are not changed. Full means that degrees of
freedom associated with dirichlet boundary conditions are set to zero in the vector. Reduced means that degrees of
freedom associated with dirichlet boundary conditions are removed and the returned vector has reduced size.
* \see const Eigen::VectorXd& getVector(const FERequirement& feRequirements, VectorAffordance affordance,
EnforcingDBCOption qt)
* \return Const reference to the calculated vectorial quantity.
*/
  const Eigen::VectorXd& vector() { return vector(this->enforcingDBCOption()); }

private:
  void assembleRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance,
                             Eigen::VectorXd& assemblyVec);
  Eigen::VectorXd& getRawVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);
  Eigen::VectorXd& getVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);

  Eigen::VectorXd& getReducedVectorImpl(const FERequirement& feRequirements, VectorAffordance affordance);

  Eigen::VectorXd vecRaw_{}; ///< Raw vector without changes for dirichlet degrees of freedom
  Eigen::VectorXd vec_{};    ///< Vector quantity.
  Eigen::VectorXd vecRed_{}; ///< Reduced vector quantity.
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
class SparseFlatAssembler : public VectorFlatAssembler<FEC, DV>
{
public:
  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using Base           = VectorFlatAssembler<FEC, DV>;

  using typename Base::Basis;
  using typename Base::DirichletValuesType;
  using typename Base::FEContainer;
  using typename Base::FERequirement;
  using typename Base::GlobalIndex;

  /**
   * \brief Constructor for SparseFlatAssembler.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  SparseFlatAssembler(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : VectorFlatAssembler<FEContainer, DirichletValuesType>(std::forward<FEContainer>(fes), dirichletValues) {}

  using GridView = typename Basis::GridView; ///< Type of the grid view.

  /**
   * \brief Calculates the matrix quantity requested by feRequirements and the affordance.
   * For EnforcingDBCOption::Full a zero is written on fixed degrees of freedom rows and columns, and a one is written
   * on the diagonal. For EnforcingDBCOption::Raw the untouched matrix is returned.
   * For EnforcingDBCOption::Reduced the matrix is reduced in size by removing the fixed degrees of freedom.

    \param feRequirements Reference to the finite element requirements.
   * \param affordance The matrix affordance
   * \param qt The EnforcingDBCOption
   * \return Const reference to the modified sparse matrix quantity.
   */
  const Eigen::SparseMatrix<double>& matrix(const FERequirement& feRequirements, MatrixAffordance affordance,
                                            EnforcingDBCOption qt = EnforcingDBCOption::Full) {
    if (qt == EnforcingDBCOption::Raw) {
      return getRawMatrixImpl(feRequirements, affordance);
    } else if (qt == EnforcingDBCOption::Reduced) {
      return getReducedMatrixImpl(feRequirements, affordance);
    } else if (qt == EnforcingDBCOption::Full) {
      return getMatrixImpl(feRequirements, affordance);
    }
    __builtin_unreachable();
  }

  /**
   * \brief Calculates the matrix quantity requested by the bound feRequirements and the affordance.
   * \see const Eigen::SparseMatrix<double>& getMatrix(const FERequirement& feRequirements,MatrixAffordance affordance,
   EnforcingDBCOption qt)

   * \param qt The EnforcingDBCOption
   * \return Const reference to the modified sparse matrix quantity.
   */
  const Eigen::SparseMatrix<double>& matrix(EnforcingDBCOption qt) {
    return matrix(this->requirement(), this->affordanceCollection().matrixAffordance(), qt);
  }

  /**
 * \brief Calculates the matrix quantity requested by the bound feRequirements, the affordance and the
enforcingDBCOption.
 * \see const Eigen::SparseMatrix<double>& getMatrix(const FERequirement& feRequirements,MatrixAffordance affordance,
 EnforcingDBCOption qt)

 * \return Const reference to the modified sparse matrix quantity.
 */
  const Eigen::SparseMatrix<double>& matrix() { return matrix(this->enforcingDBCOption()); }

private:
  void assembleRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance,
                             Eigen::SparseMatrix<double>& assemblyMat);
  Eigen::SparseMatrix<double>& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  Eigen::SparseMatrix<double>& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  Eigen::SparseMatrix<double>& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);

  /** Calculates the non-zero entries in the full sparse matrix and passes them to the underlying Eigen sparse matrix.
   */
  void createOccupationPattern(Eigen::SparseMatrix<double>& assemblyMat);

  /** Calculates the non-zero entries in the sparse matrix and passes them to the underlying Eigen sparse matrix.
   * The size of the matrix has the size of the free degrees of freedom. */
  void createReducedOccupationPattern(Eigen::SparseMatrix<double>& assemblyMat);

  /** Saves the degree of freedom indices of each element in the vector elementLinearIndices. */
  void createLinearDOFsPerElement(Eigen::SparseMatrix<double>& assemblyMat);

  /** Saves the degree of freedom indices of each element in the vector elementLinearIndices but excludes fixed
   * degrees of freedom. */
  void createLinearDOFsPerElementReduced(Eigen::SparseMatrix<double>& assemblyMat);

  /** Pre-processes the raw sparse matrix before assembly. */
  void preProcessSparseMatrix(Eigen::SparseMatrix<double>& assemblyMat);

  /** Pre-processes the reduced sparse matrix before assembly. */
  void preProcessSparseMatrixReduced(Eigen::SparseMatrix<double>& assemblyMat);

  Eigen::SparseMatrix<double> spMatRaw_;     ///< Raw sparse matrix without changes for dirichlet degrees of freedom.
  Eigen::SparseMatrix<double> spMat_;        ///< Sparse matrix.
  Eigen::SparseMatrix<double> spMatReduced_; ///< Reduced sparse matrix.
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
class DenseFlatAssembler : public VectorFlatAssembler<FEC, DV>
{
public:
  using FEContainerRaw = std::remove_cvref_t<FEC>;     ///< Type of the raw finite element container.
  using Base           = VectorFlatAssembler<FEC, DV>; ///< Type alias for the base class.

  using typename Base::Basis;               ///< Type of the basis.
  using typename Base::DirichletValuesType; ///< Type of the Dirichlet values.
  using typename Base::FEContainer;         ///< Type of the finite element container.
  using typename Base::FERequirement;       ///< Type of the finite element requirement.
  using typename Base::GlobalIndex;         ///< Type of the global index.

  /**
   * \brief Constructor for DenseFlatAssembler.
   *
   * \param fes Finite element container.
   * \param dirichletValues Reference to Dirichlet values.
   */
  explicit DenseFlatAssembler(FEContainer&& fes, const DirichletValuesType& dirichletValues)
      : VectorFlatAssembler<FEContainer, DirichletValuesType>(std::forward<FEContainer>(fes), dirichletValues) {}

  /**
   * \brief  Calculates the matrix quantity requested by feRequirements and the affordance.
   * For EnforcingDBCOption::Full a zero is written on fixed degrees of freedom rows and columns, and a one is written
   * on the diagonal. For EnforcingDBCOption::Raw the untouched matrix is returned.
   * For EnforcingDBCOption::Reduced the matrix is reduced in size by removing the fixed degrees of freedom.
   *
   * \param qt The EnforcingDBCOption

   * \return Reference to the raw dense matrix quantity.
   */
  const Eigen::MatrixXd& matrix(const FERequirement& feRequirements, MatrixAffordance affordance,
                                EnforcingDBCOption qt = EnforcingDBCOption::Full) {
    if (qt == EnforcingDBCOption::Raw) {
      return getRawMatrixImpl(feRequirements, affordance);
    } else if (qt == EnforcingDBCOption::Reduced) {
      return getReducedMatrixImpl(feRequirements, affordance);
    } else if (qt == EnforcingDBCOption::Full) {
      return getMatrixImpl(feRequirements, affordance);
    }
    __builtin_unreachable();
  }

  /**
* \brief  Calculates the matrix quantity requested by the bound  feRequirements and the affordance.
* For EnforcingDBCOption::Full a zero is written on fixed degrees of freedom rows and columns, and a one is written
* on the diagonal. For EnforcingDBCOption::Raw the untouched matrix is returned.
* For EnforcingDBCOption::Reduced the matrix is reduced in size by removing the fixed degrees of freedom.
*
* \param feRequirements Reference to the finite element requirements.
* \param affordance The matrix affordance
* \param qt The EnforcingDBCOption

* \return Reference to the raw dense matrix quantity.
*/
  const Eigen::MatrixXd& matrix(EnforcingDBCOption qt) {
    return matrix(this->requirement(), this->affordanceCollection().matrixAffordance(), qt);
  }

  /**
   * \brief  Calculates the matrix quantity requested by the bound  feRequirements, the affordance and the
enforcingDBCOption.
   * \see const Eigen::MatrixXd& matrix(EnforcingDBCOption qt)
   * \return Reference to the dense matrix quantity.
   */
  const Eigen::MatrixXd& matrix() { return matrix(this->enforcingDBCOption()); }

private:
  void assembleRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance,
                             Eigen::MatrixXd& assemblyMat);
  Eigen::MatrixXd& getRawMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  Eigen::MatrixXd& getMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);
  Eigen::MatrixXd& getReducedMatrixImpl(const FERequirement& feRequirements, MatrixAffordance affordance);

  Eigen::MatrixXd matRaw_{}; ///< Raw dense matrix for assembly.
  Eigen::MatrixXd mat_{};    ///< Dense matrix quantity.
  Eigen::MatrixXd matRed_{}; ///< Reduced dense matrix quantity.
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
