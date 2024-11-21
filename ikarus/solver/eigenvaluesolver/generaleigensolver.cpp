// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "generaleigensolver.hh"

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Ikarus {

template class GeneralSymEigenSolver<EigenSolverTypeTag::Spectra, Eigen::SparseMatrix<double>>;
template class GeneralSymEigenSolver<EigenSolverTypeTag::Spectra, Eigen::MatrixX<double>>;
template class GeneralSymEigenSolver<EigenSolverTypeTag::Eigen, Eigen::MatrixX<double>>;

template class PartialGeneralSymEigenSolver<Eigen::MatrixX<double>>;

} // namespace Ikarus
