//
// Created by lex on 04/02/2022.
//

#pragma once
#include <ranges>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

#include <Eigen/Core>

namespace Ikarus {
  template <typename DuneLocalBasis>
  class LocalBasis {
  public:
    LocalBasis(const DuneLocalBasis& p_basis) : duneLocalBasis{&p_basis} {}
    LocalBasis() = default;

    static constexpr int gridDim = DuneLocalBasis::Traits::dimDomain;
    using DomainType             = typename DuneLocalBasis::Traits::DomainType;
    using RangeType              = typename DuneLocalBasis::Traits::RangeType;
    using DomainFieldType        = typename DuneLocalBasis::Traits::DomainFieldType;
    using RangeFieldType         = typename DuneLocalBasis::Traits::RangeFieldType;
    using JacobianType           = typename DuneLocalBasis::Traits::JacobianType;

    template <typename Derived>
    void evaluateFunction(const DomainType& local, Eigen::PlainObjectBase<Derived>& N) const {
      duneLocalBasis->evaluateFunction(local, Ndune);
      N.resize(Ndune.size(), 1);
      N.setZero();
      for (size_t i = 0; i < Ndune.size(); ++i)
        N[i] = Ndune[i][0];
    }

    template <typename Derived>
    void evaluateJacobian(const DomainType& local, Eigen::PlainObjectBase<Derived>& dN) const {
      duneLocalBasis->evaluateJacobian(local, dNdune);
      dN.setZero();
      dN.resize(dNdune.size(), gridDim);

      for (auto i = 0U; i < dNdune.size(); ++i)
        for (int j = 0; j < gridDim; ++j)
          dN(i, j) = dNdune[i][0][j];
    }

    template <typename Derived1, typename Derived2>
    void evaluateFunctionAndJacobian(const DomainType& local, Eigen::PlainObjectBase<Derived1>& N,
                                     Eigen::PlainObjectBase<Derived2>& dN) const {
      evaluateFunction(local, N);
      evaluateJacobian(local, dN);
    }

    unsigned int size() { return duneLocalBasis->size(); }

    template <typename IntegrationRule, typename... Ints>
    requires std::conjunction_v<std::is_convertible<int, Ints>...>
    void bind(IntegrationRule&& p_rule, Ints&&... ints) {
      rule             = p_rule;
      boundDerivatives = std::vector<int>({std::forward<Ints>(ints)...});
      Nbound           = std::make_optional(std::vector<Eigen::VectorX<RangeFieldType>>{});
      dNbound          = std::make_optional(std::vector<Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>>{});
      dNbound.value().resize(rule.value().size());
      Nbound.value().resize(rule.value().size());
      std::ranges::sort(boundDerivatives.value());

      for (int i = 0; auto& gp : rule.value()) {
        Eigen::Matrix<double, Eigen::Dynamic, gridDim> dN;
        if (std::ranges::binary_search(boundDerivatives.value(), 0)) evaluateFunction(gp.position(), Nbound.value()[i]);
        if (std::ranges::binary_search(boundDerivatives.value(), 1))
          evaluateJacobian(gp.position(), dNbound.value()[i]);
        ++i;
      }
    }

    const auto& getFunction(long unsigned i) const {
      if (not Nbound) throw std::logic_error("You have to bind the basis first");
      return Nbound.value()[i];
    }

    const auto& getJacobian(long unsigned i) const {
      if (not dNbound) throw std::logic_error("You have to bind the basis first");
      return dNbound.value()[i];
    }

    bool isBound()const {return (dNbound and Nbound);}

    struct FunctionAndJacobian {
      long unsigned index;
      const Dune::QuadraturePoint<DomainFieldType, gridDim>& ip;
      const Eigen::VectorX<RangeFieldType>& N;
      const Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& dN;
    };
    auto viewOverFunctionAndJacobian() const {
      assert(Nbound.value().size() == dNbound.value().size()
             && "Number of intergrationpoint evaluations does not match.");
      if (Nbound and dNbound)
        return std::views::iota(0UL, Nbound.value().size()) | std::views::transform([&](auto&& i_) {
                 return FunctionAndJacobian(i_, rule.value()[i_], getFunction(i_), getJacobian(i_));
               });
      else {
        assert(false && "You need to call bind first");
        __builtin_unreachable();
      }
    }

  private:
    mutable std::vector<JacobianType> dNdune;
    mutable std::vector<RangeType> Ndune;
    DuneLocalBasis const* duneLocalBasis;
    std::optional<std::vector<int>> boundDerivatives;
    std::optional<std::vector<Eigen::VectorX<RangeFieldType>>> Nbound;
    std::optional<std::vector<Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>>> dNbound;
    std::optional<Dune::QuadratureRule<DomainFieldType, gridDim>> rule;
  };

}  // namespace Ikarus