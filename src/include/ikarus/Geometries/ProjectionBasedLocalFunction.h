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
#include "ikarus/utils/namedType.h"
#include <ikarus/LocalBasis/localBasis.h>
#include <ikarus/utils/LinearAlgebraHelper.h>

namespace Ikarus {
  namespace DerivativeDirections {
    static struct Coeffs {
    } coeffs;
    static struct Spatial { } spatial; }  // namespace DerivativeDirections

  template <typename... Args>
  struct CoeffIndices {
    std::array<std::common_type_t<std::remove_cvref_t<Args>...>, sizeof...(Args)> args;
  };

  template <typename... Args>
  auto coeffIndices(Args&&... args) {
    return CoeffIndices<Args&&...>({std::forward<Args>(args)...});
  }

  template <typename... Args>
  struct Wrt {
    std::tuple<std::remove_cvref_t<Args>...> args;
  };

  template <typename... Args>
  auto wrt(Args&&... args) {
    return Wrt<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
  }

  template <typename... Args>
  struct Along {
    std::tuple<Args...> args;
  };

  template <typename... Args>
  auto along(Args&&... args) {
    return Along<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
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
          if constexpr (std::is_same_v<DerivativeDirections::Coeffs, CurrentType>)
            ++counter.coeffDerivatives;
          else if constexpr (std::is_same_v<DerivativeDirections::Spatial, CurrentType>)
            ++counter.spatialDerivatives;
          else
            static_assert(std::is_same_v<DerivativeDirections::Coeffs,
                                         CurrentType> or std::is_same_v<DerivativeDirections::Spatial, CurrentType>,
                          "You are only allowed to pass coeffs and spatial struct");
        });
    return counter;
  }

  template <typename DuneBasis, typename CoeffContainer>
  class ProjectionBasedLocalFunction {
  public:
    using DomainType = typename DuneBasis::Traits::DomainType;
    ProjectionBasedLocalFunction(const Ikarus::LocalBasis<DuneBasis>& basis_, const CoeffContainer& coeffs_)
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
    using FieldMat           = Eigen::Matrix<ctype, coeffdimension, coeffdimension>;

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

    template <typename... Args, typename... Indices>
    requires((countInTuple<Wrt<Args...>>().coeffDerivatives == 1)
             and (countInTuple<Wrt<Args...>>().coeffDerivatives
                  == 1)) auto evaluateDerivative(long unsigned gpIndex, const Manifold& val, Wrt<Args...>&& args,
                                                 CoeffIndices<Indices...>&& coeffsIndices) const {
      constexpr Counter counter = countInTuple<Wrt<Args...>>();
      if constexpr (counter.coeffDerivatives == 0 and counter.spatialDerivatives == 1) {
        return evaluateDerivative(gpIndex, val);
      } else if constexpr (counter.coeffDerivatives == 1 and counter.spatialDerivatives == 0) {
        return evaluateDerivativeWRTCoeffs(gpIndex, val, coeffsIndices.args[0]);
      } else if constexpr (counter.coeffDerivatives == 1 and counter.spatialDerivatives == 1) {
        return evaluateDerivativeWRTCoeffsANDSpatial(gpIndex, val, coeffsIndices.args[0]);
      }
    }

    template <typename... Args, typename... AlongArgs, typename... Indices>
    requires((countInTuple<Wrt<Args...>>().coeffDerivatives == 2)
             and (countInTuple<Wrt<Args...>>().coeffDerivatives
                  == 2)) auto evaluateDerivative(long unsigned gpIndex, const Manifold& val, Wrt<Args...>&& args,
                                                 Along<AlongArgs...>&& along,
                                                 CoeffIndices<Indices...>&& coeffsIndices) const {
      constexpr Counter counter = countInTuple<Wrt<Args...>>();
      if constexpr (counter.coeffDerivatives == 2 and counter.spatialDerivatives == 0) {
        return evaluateSecondDerivativeWRTCoeffs(gpIndex, val, std::get<0>(along.args), coeffsIndices.args);
      }else if constexpr (counter.coeffDerivatives == 2 and counter.spatialDerivatives == 1) {
        return evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatial(gpIndex, val, std::get<0>(along.args), coeffsIndices.args);
      }
    }

    auto evaluateDerivativeWRTCoeffsANDSpatial(const long unsigned gpIndex, const Manifold& val,
                                               int coeffsIndex) const {
      const GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(gpIndex));
      const Jacobian J   = evaluateEmbeddingJacobianImpl(gpIndex);
      const auto Pm      = Manifold::derivativeOfProjectionWRTposition(valE);
      std::array<FieldMat, gridDim> Warray;
      for (int i = 0; i < gridDim; ++i) {
        const auto Qi = Manifold::secondDerivativeOfProjectionWRTposition(valE, J.col(i));
        Warray[i]     = Qi * basis.getFunction(gpIndex)[coeffsIndex] + Pm * basis.getJacobian(gpIndex)(coeffsIndex, i);
      }

      return Warray;
    }

    auto evaluateDerivativeWRTCoeffs(const long unsigned gpIndex, const Manifold& val, int coeffsIndex) const {
      GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(gpIndex));
      return (Manifold::derivativeOfProjectionWRTposition(valE) * basis.getFunction(gpIndex)[coeffsIndex]).eval();
    }

    auto evaluateSecondDerivativeWRTCoeffs(const long unsigned gpIndex, const Manifold& val,
                                           const Eigen::Vector<ctype, coeffdimension>& along,
                                           const std::array<size_t, gridDim>& coeffsIndex) const {
      const GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(gpIndex));
      FieldMat Snn       = Manifold::secondDerivativeOfProjectionWRTposition(valE, along)
                     * basis.getFunction(gpIndex)[coeffsIndex[0]] * basis.getFunction(gpIndex)[coeffsIndex[1]];

      return Snn;
    }

    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatial(const long unsigned gpIndex, const Manifold& val,
                                                            const Eigen::Vector<ctype, coeffdimension>& along,
                                                            const std::array<size_t, gridDim>& coeffsIndex) const {
      const GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(gpIndex));
      const Jacobian J   = evaluateEmbeddingJacobianImpl(gpIndex);
      const auto S       = Manifold::secondDerivativeOfProjectionWRTposition(valE, along);
      std::array<FieldMat, gridDim> ChiArray;
      for (int i = 0; i < gridDim; ++i) {
        const auto chi    = Manifold::thirdDerivativeOfProjectionWRTposition(valE, along, J.col(i));
        const auto& NI    = basis.getFunction(gpIndex)[coeffsIndex[0]];
        const auto& NJ    = basis.getFunction(gpIndex)[coeffsIndex[1]];
        const auto& dNIdi = basis.getJacobian(gpIndex)(coeffsIndex[0], i);
        const auto& dNJdi = basis.getJacobian(gpIndex)(coeffsIndex[1], i);
        ChiArray[i]       = chi * NI * NJ + S * (dNIdi * NJ + dNJdi * NI);
      }

      return ChiArray;
    }

    template <typename... Args>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives == 0) Jacobian
        evaluateDerivative(long unsigned i, const Manifold& val, Wrt<Args...>&& args)
    const {
      constexpr Counter counter = countInTuple<Wrt<Args...>>();
      static_assert(counter.spatialDerivatives < 2, "This currently only supports first order spatial derivatives");
      if constexpr (counter.spatialDerivatives == 1) {
        return evaluateFirstSpatialDerivative(i, val);
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

    template <typename... Args, typename... Indices>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives
             == 1) auto evaluateDerivative(long unsigned gpIndex, Wrt<Args...>&& args,
                                           CoeffIndices<Indices...>&& coeffsIndices) const {
      Manifold val = evaluateFunctionImpl(basis.getFunction(gpIndex));
      return evaluateDerivative(gpIndex, val, std::forward<Wrt<Args...>>(args),
                                std::forward<CoeffIndices<Indices...>>(coeffsIndices));
    }

    template <typename... Args, typename... AlongArgs, typename... Indices>
    requires(countInTuple<Wrt<Args...>>().coeffDerivatives
             == 2) auto evaluateDerivative(long unsigned gpIndex, Wrt<Args...>&& args, Along<AlongArgs...>&& along,
                                           CoeffIndices<Indices...>&& coeffsIndices) const {
      Manifold val = evaluateFunctionImpl(basis.getFunction(gpIndex));
      return evaluateDerivative(gpIndex, val, std::forward<Wrt<Args...>>(args),
                                std::forward<Along<AlongArgs...>>(along),
                                std::forward<CoeffIndices<Indices...>>(coeffsIndices));
    }

    auto viewOverIntegrationPoints()
    {
      return basis.viewOverIntegrationPoints();
    }
  private:
    Jacobian evaluateFirstSpatialDerivative(long unsigned i, const Manifold& val) const {
      assert(basis.isBound() && "You have to bind the basis first");
      Jacobian J = evaluateEmbeddingJacobianImpl(i);

      GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(i));
      return Manifold::derivativeOfProjectionWRTposition(valE) * J;
    }

    Jacobian evaluateFirstSpatialDerivative(long unsigned i) const {
      Manifold val = evaluateFunctionImpl(basis.getFunction(i));
      return evaluateFirstSpatialDerivative(i, val);
    }
    Manifold evaluateFunctionImpl(const Eigen::VectorXd& N) const { return Manifold(evaluateEmbeddingFunctionImpl(N)); }
    Jacobian evaluateEmbeddingJacobianImpl(long unsigned i) const {
      assert(basis.isBound() && "You have to bind the basis first");
      const auto& dN = basis.getJacobian(i);
      Jacobian J     = coeffsAsMat * dN;
      return J;
    }

    GlobalE evaluateEmbeddingFunctionImpl(const Eigen::VectorXd& N) const { return coeffsAsMat * N; }
    const Ikarus::LocalBasis<DuneBasis>& basis;
    const CoeffContainer& coeffs;
    const decltype(Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

}  // namespace Ikarus
