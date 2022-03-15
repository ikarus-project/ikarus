//
// Created by Alex on 21.04.2021.
//

#pragma once

#include <concepts>
#include <iostream>
#include <span>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Interpolators/Interpolator.h"
#include <ikarus/LocalBasis/localBasis.h>
#include <ikarus/utils/LinearAlgebraHelper.h>

namespace Ikarus {

  struct coeffs {};
  struct spatial {};

  template <typename... Args>
  struct Wrt {
    std::tuple<std::remove_cvref_t<Args>...> args;
  };

  template <typename... Args>
  auto wrt(Args&&... args) {
    return Wrt<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
  }


  struct Counter {
    int coeffDerivatives;
    int spatialDerivatives;
  };
  template <typename WrtType>
  consteval Counter countInTuple() {
    Counter counter{0, 0};
    Dune::Hybrid::forEach(
        Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(WrtType::args)>>()), [&](auto i) {
          using CurrentType = std::tuple_element_t<i, decltype(WrtType::args)>;
          if constexpr (std::is_same_v<coeffs, CurrentType>)
            ++counter.coeffDerivatives;
          else if constexpr (std::is_same_v<spatial, CurrentType>)
            ++counter.spatialDerivatives;
          else
            static_assert(std::is_same_v<coeffs, CurrentType> or std::is_same_v<spatial, CurrentType>,
                          "You are only allowed to pass coeffs and spatial struct");
        });
    return counter;
  }

  template <typename DuneBasis, typename CoeffContainer>
  class ProjectionBasedLocalFunction {
  public:
    using DomainType = typename DuneBasis::Traits::DomainType;
    ProjectionBasedLocalFunction(Ikarus::LocalBasis<DuneBasis>& basis_, const CoeffContainer& coeffs_)
        : basis{basis_}, coeffs{coeffs_}, coeffsAsMat{Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)} {}
    /** \brief Type used for coordinates */
    using ctype = typename CoeffContainer::value_type::ctype;
    /** \brief Dimension of the coeffs */
    static constexpr int coeffdimension = CoeffContainer::value_type::valueSize;

    /** \brief Dimension of the grid */
    static constexpr int gridDim = Ikarus::LocalBasis<DuneBasis>::gridDim;

    /** \brief Type for coordinate vector in world space */
    using Manifold = typename CoeffContainer::value_type;
    using GlobalE  = decltype(std::declval<typename CoeffContainer::value_type>().getValue());

    /** \brief Type for the transposed Jacobian matrix */
    using JacobianTransposed = Eigen::Matrix<ctype, gridDim, coeffdimension>;
    using Jacobian           = Eigen::Matrix<ctype, coeffdimension, gridDim>;

    /** \brief Type for the transposed inverse Jacobian matrix */
    using JacobianInverseTransposed = Eigen::Matrix<ctype, coeffdimension, gridDim>;
    using JacobianInverse           = Eigen::Matrix<ctype, gridDim, coeffdimension>;

    Manifold evaluateFunction(long unsigned i) {
      assert(basis.isBound() && "You have to bind the basis first");
      const auto& N = basis.getFunction(i);
      return evaluateFunctionImpl(N);
    }

    Manifold evaluateFunction(const DomainType& local) {
      Eigen::VectorXd N;
      basis.evaluateFunction(local, N);
      return evaluateFunctionImpl(N);
    }

    Jacobian evaluateDerivative(long unsigned i, const Manifold& val) const {
      assert(basis.isBound() && "You have to bind the basis first");
      const auto& dN = basis.getJacobian(i);
      Jacobian J     = coeffsAsMat * dN;

      GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(i));
      return Manifold::derivativeOfProjection(valE) * J;
    }

    Jacobian evaluateDerivative(long unsigned i) const {
      Manifold val = evaluateFunctionImpl(basis.getFunction(i));
      return evaluateDerivative(i, val);
    }

    template <typename... Args>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives > 0) auto evaluateDerivative(
        long unsigned gpIndex, const Manifold& val, Wrt<Args...>&& args,
        std::array<size_t, countInTuple<Wrt<Args...>>().coeffDerivatives>&& coeffsIndices) const {
      constexpr Counter counter = countInTuple<Wrt<Args...>>();

      if constexpr (counter.coeffDerivatives == 0 and counter.spatialDerivatives == 1) {
        return evaluateDerivative(gpIndex, val);
      } else if constexpr (counter.coeffDerivatives == 1 and counter.spatialDerivatives == 0) {
        GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(gpIndex));
        return (Manifold::derivativeOfProjection(valE) * basis.getFunction(gpIndex)[coeffsIndices[0]]).eval();
      }else if constexpr (counter.coeffDerivatives == 1 and counter.spatialDerivatives == 1) {
        GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(gpIndex)); //impl W
        return (Manifold::derivativeOfProjection(valE) * basis.getFunction(gpIndex)[coeffsIndices[0]]).eval();
      }
    }

    template <typename... Args>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives == 0) Jacobian
        evaluateDerivative(long unsigned i, const Manifold& val, Wrt<Args...>&& args)
    const {
      constexpr Counter counter = countInTuple<Wrt<Args...>>();
      static_assert(counter.spatialDerivatives < 2, "This currently only supports first order spatial derivatives");
      if constexpr (counter.spatialDerivatives == 1) {
        return evaluateDerivative(i, val);
      } else
        __builtin_unreachable();
    }

    template <typename... Args>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives == 0) Jacobian
        evaluateDerivative(long unsigned gpIndex, Wrt<Args...>&& args)
    const {
      Manifold val = evaluateFunctionImpl(basis.getFunction(gpIndex));
      return evaluateDerivative(gpIndex, val, std::forward<Wrt<Args...>>(args));
    }

    template <typename... Args>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives > 0) Jacobian
        evaluateDerivative(long unsigned gpIndex, Wrt<Args...>&& args,
                           std::array<size_t, countInTuple<Wrt<Args...>>().coeffDerivatives>&& coeffsIndices)
    const {
      Manifold val = evaluateFunctionImpl(basis.getFunction(gpIndex));
      return evaluateDerivative(
          gpIndex, val, std::forward<Wrt<Args...>>(args),
          std::forward<std::array<size_t, countInTuple<Wrt<Args...>>().coeffDerivatives>>(coeffsIndices));
    }

  private:
    Manifold evaluateFunctionImpl(const Eigen::VectorXd& N) const { return Manifold(evaluateEmbeddingFunctionImpl(N)); }

    GlobalE evaluateEmbeddingFunctionImpl(const Eigen::VectorXd& N) const { return coeffsAsMat * N; }
    const Ikarus::LocalBasis<DuneBasis>& basis;
    const CoeffContainer& coeffs;
    const decltype(Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

}  // namespace Ikarus
