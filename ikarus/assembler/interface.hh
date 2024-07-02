// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Defines the interface for scalar, vector and matrix assemblers
 */

#pragma once

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/dirichletvalues.hh>

namespace Ikarus {

template <typename SA, typename FEC, typename DV, typename ST>
class ScalarAssembler
{
public:
  using ScalarAssemblerType = SA;
  using ScalarType          = ST;
  using DirichletValuesType = DV;
  using FEContainer         = FEC;
  using FEContainerRaw      = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using FERequirement       = typename FEContainerRaw::value_type::Requirement;

  /**
   * \brief Calculates the scalar quantity requested by feRequirements and affordance.
   *
   * \param feRequirements Reference to the finite element requirements.
   * \param affordance The scalar affordance
   * \return Const reference to the calculated scalar quantity.
   */
  const ScalarType& scalar(const FERequirement& feRequirements, ScalarAffordance affordance) {
    return underlying().getScalarImpl(feRequirements, affordance);
  }

  /**
   * \brief Calculates the scalar quantity requested by the bound feRequirements and returns a reference.
   *
   * \return Const reference to the calculated scalar quantity.
   */
  const ScalarType& scalar() {
    return underlying().getScalarImpl(underlying().requirement(),
                                      underlying().affordanceCollection().scalarAffordance());
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const ScalarAssemblerType&>(*this); }
  auto& underlying() { return static_cast<ScalarAssemblerType&>(*this); }
};

template <typename VA, typename FEC, typename DV, typename VT>
class VectorAssembler
{
public:
  using VectorAssemblerType = VA;
  using VectorType          = VT;

  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using FERequirement  = typename FEContainerRaw::value_type::Requirement;
  ///< Type of the finite element requirement.
  using GlobalIndex = typename FEContainerRaw::value_type::GlobalIndex; ///< Type of the global index.

  using DirichletValuesType = DV;
  using FEContainer         = FEC;

  /**
     * \brief Calculates the vectorial quantity requested by the  feRequirements and the affordance.
     Depending on the requested DBCOption, the raw, reduced or full vector is returned.
      Raw means the degrees of freedom associated with dirichlet boundary conditions are not changed.
      Full means that degrees of freedom associated with dirichlet boundary conditions are set to zero in the vector.
      Reduced means that degrees of freedom associated with dirichlet boundary conditions are removed and the returned
     vector has reduced size.
     *
     * \param feRequirements Reference to the finite element requirements.
     * \param affordance The vector affordance
     * \param dbcOption The DBCOption
     * \return Const reference to the calculated vectorial quantity.
     */
  const VectorType& vector(const FERequirement& feRequirements, VectorAffordance affordance,
                           DBCOption dbcOption = DBCOption::Full) {
    if (dbcOption == DBCOption::Raw) {
      return underlying().getRawVectorImpl(feRequirements, affordance);
    } else if (dbcOption == DBCOption::Reduced) {
      return underlying().getReducedVectorImpl(feRequirements, affordance);
    } else if (dbcOption == DBCOption::Full) {
      return underlying().getVectorImpl(feRequirements, affordance);
    }
    __builtin_unreachable();
  }

  /**
 * \brief Calculates the vectorial quantity requested by the bound feRequirements and the affordance.
 Depending on the requested DBCOption, the raw, reduced or full vector is returned.
  Raw means the degrees of freedom associated with dirichlet boundary conditions are not changed.
  Full means that degrees of freedom associated with dirichlet boundary conditions are set to zero in the vector.
  Reduced means that degrees of freedom associated with dirichlet boundary conditions are removed and the returned
 vector has reduced size.
 * \param dbcOption The DBCOption
 * \return Const reference to the calculated vectorial quantity.
 */
  const VectorType& vector(DBCOption dbcOption) {
    return vector(underlying().requirement(), underlying().affordanceCollection().vectorAffordance(), dbcOption);
  }

  /**
* \brief Calculates the vectorial quantity requested by the bound feRequirements,  the affordance and the
dBCOption. Depending on the DBCOption, the raw, reduced or full vector is returned. Raw
means the degrees of freedom associated with dirichlet boundary conditions are not changed. Full means that degrees of
freedom associated with dirichlet boundary conditions are set to zero in the vector. Reduced means that degrees of
freedom associated with dirichlet boundary conditions are removed and the returned vector has reduced size.
* \return Const reference to the calculated vectorial quantity.
*/
  const VectorType& vector() { return vector(underlying().dBCOption()); }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const VectorAssemblerType&>(*this); }
  auto& underlying() { return static_cast<VectorAssemblerType&>(*this); }
};

template <typename MA, typename FEC, typename DV, typename MT>
class MatrixAssembler
{
public:
  using MatrixAssemblerType = MA;
  using MatrixType          = MT;

  using FEContainerRaw = std::remove_cvref_t<FEC>; ///< Type of the raw finite element container.
  using FERequirement  = typename FEContainerRaw::value_type::Requirement;
  ///< Type of the finite element requirement.
  using GlobalIndex = typename FEContainerRaw::value_type::GlobalIndex; ///< Type of the global index.

  using DirichletValuesType = DV;
  using FEContainer         = FEC;

  /**
   * \brief Calculates the matrix quantity requested by feRequirements and the affordance.
   * For DBCOption::Full a zero is written on fixed degrees of freedom rows and columns, and a one is written
   * on the diagonal. For DBCOption::Raw the untouched matrix is returned.
   * For DBCOption::Reduced the matrix is reduced in size by removing the fixed degrees of freedom.

    \param feRequirements Reference to the finite element requirements.
   * \param affordance The matrix affordance
   * \param dbcOption The DBCOption
   * \return Const reference to the modified sparse matrix quantity.
   */
  const MatrixType& matrix(const FERequirement& feRequirements, MatrixAffordance affordance,
                           DBCOption dbcOption = DBCOption::Full) {
    if (dbcOption == DBCOption::Raw) {
      return underlying().getRawMatrixImpl(feRequirements, affordance);
    } else if (dbcOption == DBCOption::Reduced) {
      return underlying().getReducedMatrixImpl(feRequirements, affordance);
    } else if (dbcOption == DBCOption::Full) {
      return underlying().getMatrixImpl(feRequirements, affordance);
    }
    __builtin_unreachable();
  }

  /**
   * \brief Calculates the matrix quantity requested by the bound feRequirements and the affordance.
   * \see const Eigen::SparseMatrix<double>& matrix(const FERequirement& feRequirements,MatrixAffordance affordance,
   DBCOption dbcOption)

   * \param dbcOption The DBCOption
   * \return Const reference to the modified sparse matrix quantity.
   */
  const MatrixType& matrix(DBCOption dbcOption) {
    return matrix(underlying().requirement(), underlying().affordanceCollection().matrixAffordance(), dbcOption);
  }

  /**
 * \brief Calculates the matrix quantity requested by the bound feRequirements, the affordance and the
dBCOption.
 * \see const Eigen::SparseMatrix<double>& matrix(const FERequirement& feRequirements,MatrixAffordance affordance,
 DBCOption dbcOption)

 * \return Const reference to the modified sparse matrix quantity.
 */
  const MatrixType& matrix() { return matrix(underlying().dBCOption()); }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const MatrixAssemblerType&>(*this); }
  auto& underlying() { return static_cast<MatrixAssemblerType&>(*this); }
};
} // namespace Ikarus
