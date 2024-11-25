// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file dynamics.hh
 * \brief Helper for
 */

#pragma once
#include <dune/vtk/pvdwriter.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics::LumpingSchemes {

struct RowSumLumping
{
  template <typename ScalarType = double>
  void operator()(const auto&, const auto&, auto, auto, Eigen::SparseMatrix<ScalarType>& mat) {
    for (auto i : Dune::range(mat.rows())) {
      auto sum           = mat.row(i).sum();
      mat.coeffRef(i, i) = sum;
    }
    // Deletes all entries expect for main diagonal (entriese are really getting deleted)
    mat.prune([](int i, int j, auto) { return i == j; });
  }

  template <typename ScalarType = double>
  void operator()(const auto&, const auto&, auto, auto, Eigen::MatrixX<ScalarType>& mat) {
    for (auto i : Dune::range(mat.rows())) {
      auto sum  = mat.row(i).sum();
      mat(i, i) = sum;
    }
    mat = mat.diagonal().asDiagonal();
  }
};

} // namespace Ikarus::Dynamics::LumpingSchemes
