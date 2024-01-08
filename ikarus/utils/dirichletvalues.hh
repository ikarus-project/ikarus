// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file dirichletvalues.hh
 * @brief Definition of DirichletValues class for handling Dirichlet boundary conditions.
 *
 */

#pragma once
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <Eigen/Core>

#include <autodiff/forward/real/real.hpp>

#include <ikarus/utils/concepts.hh>

namespace Dune {
  template <class K, int N>
  class FieldVector;
}

namespace Ikarus {

  /**
   * @brief Class for handling Dirichlet boundary conditions in Ikarus.
   * \ingroup  utils
   * \details The DirichletValues class provides functionalities for fixing degrees of freedom and storing inhomogeneous
   * Dirichlet boundary conditions. It supports fixing degrees of freedom using various callback functions and
   * stores functions for inhomogeneous Dirichlet boundary conditions.
   *
   * @tparam Basis_ Type of the finite element basis
   * @tparam FlagsType_ Type for storing Dirichlet flags (default is std::vector<bool>)
   */
  template <typename Basis_, typename FlagsType_ = std::vector<bool>>
  class DirichletValues {
  public:
    using Basis                         = std::remove_cvref_t<Basis_>;
    using FlagsType                     = FlagsType_;
    static constexpr int worldDimension = Basis::GridView::dimensionworld;
    using BackendType                   = decltype(Dune::Functions::istlVectorBackend(std::declval<FlagsType&>()));
    explicit DirichletValues(const Basis_& p_basis) : basis_{p_basis}, dirichletFlagsBackend{dirichletFlags} {
      dirichletFlagsBackend.resize(basis_);
      std::fill(dirichletFlags.begin(), dirichletFlags.end(), false);
      inhomogeneousBoundaryVectorDummy.setZero(static_cast<Eigen::Index>(basis_.size()));
    }

    /**
     * @brief Function to fix (set boolean values to true or false) degrees of freedom on the boundary.
     *
     * This function takes a callback function, which will be called with the boolean vector of fixed boundary
     * degrees of freedom and the usual arguments of `Dune::Functions::forEachBoundaryDOF`.
     *
     * @param f A callback function
     */
    template <typename F>
    void fixBoundaryDOFs(F&& f) {
      if constexpr (Concepts::IsFunctorWithArgs<F, BackendType, typename Basis::MultiIndex>) {
        auto lambda = [&](auto&& indexGlobal) { f(dirichletFlagsBackend, indexGlobal); };
        Dune::Functions::forEachBoundaryDOF(basis_, lambda);
      } else if constexpr (Concepts::IsFunctorWithArgs<F, BackendType, int, typename Basis::LocalView>) {
        auto lambda = [&](auto&& localIndex, auto&& localView) { f(dirichletFlagsBackend, localIndex, localView); };
        Dune::Functions::forEachBoundaryDOF(basis_, lambda);
      } else if constexpr (Concepts::IsFunctorWithArgs<F, BackendType, int, typename Basis::LocalView,
                                                       typename Basis::GridView::Intersection>) {
        auto lambda = [&](auto&& localIndex, auto&& localView, auto&& intersection) {
          f(dirichletFlagsBackend, localIndex, localView, intersection);
        };
        Dune::Functions::forEachBoundaryDOF(basis_, lambda);
      }
    }

    /**
     * @brief Function to fix (set boolean values to true or false) degrees of freedom.
     *
     * This function takes a callback function, which will be called with the stored function space basis and
     * the boolean vector of fixed boundary degrees of freedom.
     *
     * @param f A callback function
     */
    template <typename F>
    void fixDOFs(F&& f) {
      f(basis_, dirichletFlagsBackend);
    }

    /**
     * \brief Function to fix (set boolean values to true or false) of degrees of freedom
     *
     * \param i An index indicating the DOF number to be fixed
     */
    void fixIthDOF(typename Basis::MultiIndex i) { dirichletFlagsBackend[i] = true; }

    /* \brief Returns the local basis object */
    const auto& basis() const { return basis_; }

    /* \brief Returns a boolean values, if the give multiIndex is constrained */
    template <typename MultiIndex>
    requires(not std::integral<MultiIndex>) [[nodiscard]] bool isConstrained(const MultiIndex& multiIndex) const {
      return dirichletFlagsBackend[multiIndex];
    }

    /* \brief Returns a boolean values, if the i-th degree of freedom is constrained */
    [[nodiscard]] bool isConstrained(std::size_t i) const { return dirichletFlags[i]; }

    /* \brief Returns how many degrees of freedoms are fixed */
    auto fixedDOFsize() const { return std::ranges::count(dirichletFlags, true); }

    /* \brief Returns how many degrees of freedom there are */
    auto size() const { return dirichletFlags.size(); }

    /* \brief Returns the underlying dof flag container */
    auto& container() const { return dirichletFlags; }

    /**
     * @brief Function to insert a function for inhomogeneous Dirichlet boundary conditions.
     *
     * \details This function takes a callback function, which will be called with the current coordinate vector
     * and the scalar load factor. It creates internally the first derivative of the passed function and
     * stores them simultaneously.
     *
     * @param f A callback function
     */
    template <typename F>
    void storeInhomogeneousBoundaryCondition(F&& f) {
      auto derivativeLambda = [&](const auto& globalCoord, const double& lambda) {
        autodiff::real lambdaDual = lambda;
        lambdaDual[1]             = 1;  // Setting the derivative in lambda direction to 1
        return derivative(f(globalCoord, lambdaDual));
      };
      dirichletFunctions.push_back({f, derivativeLambda});
    }

    /**
     * @brief Function to evaluate all stored inhomogeneous Dirichlet boundary functions.
     *
     * This function evaluates all stored inhomogeneous Dirichlet boundary functions at all positions where the
     * corresponding degrees of freedom are true.
     *
     * @param xIh The vector where the interpolated result should be stored
     * @param lambda The load factor
     */
    void evaluateInhomogeneousBoundaryCondition(Eigen::VectorXd& xIh, const double& lambda) {
      inhomogeneousBoundaryVectorDummy.setZero();
      xIh.resizeLike(inhomogeneousBoundaryVectorDummy);
      xIh.setZero();
      for (auto& f : dirichletFunctions) {
        interpolate(
            basis_, inhomogeneousBoundaryVectorDummy,
            [&](const auto& globalCoord) { return f.value(globalCoord, lambda); }, dirichletFlagsBackend);
        xIh += inhomogeneousBoundaryVectorDummy;
      }
    }

    /**
     * @brief Function to evaluate all stored inhomogeneous Dirichlet boundary derivative functions.
     *
     * This function evaluates all stored inhomogeneous Dirichlet boundary derivative functions at all positions
     * where the corresponding degrees of freedom are true.
     *
     * @param xIh The vector where the interpolated result should be stored
     * @param lambda The load factor
     */
    void evaluateInhomogeneousBoundaryConditionDerivative(Eigen::VectorXd& xIh, const double& lambda) {
      inhomogeneousBoundaryVectorDummy.setZero();
      xIh.resizeLike(inhomogeneousBoundaryVectorDummy);
      xIh.setZero();
      for (auto& f : dirichletFunctions) {
        interpolate(
            basis_, inhomogeneousBoundaryVectorDummy,
            [&](const auto& globalCoord) { return f.derivative(globalCoord, lambda); }, dirichletFlagsBackend);
        xIh += inhomogeneousBoundaryVectorDummy;
      }
    }

  private:
    Eigen::VectorXd inhomogeneousBoundaryVectorDummy;
    Basis basis_;
    FlagsType dirichletFlags;
    BackendType dirichletFlagsBackend;
    struct DirichletFunctions {
      using Signature = std::function<Eigen::Vector<double, worldDimension>(
          const Dune::FieldVector<double, worldDimension>&, const double&)>;
      Signature value;
      Signature derivative;
    };
    std::vector<DirichletFunctions> dirichletFunctions;
  };

}  // namespace Ikarus
