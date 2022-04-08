

namespace Ikarus {

  template <Concepts::DuneLocalBasis DuneLocalBasis>
  template <typename Derived>
  void LocalBasis<DuneLocalBasis>::evaluateFunction(const DomainType& local, Eigen::PlainObjectBase<Derived>& N) const {
    duneLocalBasis->evaluateFunction(local, Ndune);
    N.resize(Ndune.size(), 1);
    N.setZero();
    for (size_t i = 0; i < Ndune.size(); ++i)
      N[i] = Ndune[i][0];
  }

  template <Concepts::DuneLocalBasis DuneLocalBasis>
  template <typename Derived>
  void LocalBasis<DuneLocalBasis>::evaluateJacobian(const DomainType& local, Eigen::PlainObjectBase<Derived>& dN) const {
    duneLocalBasis->evaluateJacobian(local, dNdune);
    dN.setZero();
    dN.resize(dNdune.size(), gridDim);

    for (auto i = 0U; i < dNdune.size(); ++i)
      for (int j = 0; j < gridDim; ++j)
        dN(i, j) = dNdune[i][0][j];
  }

  template <Concepts::DuneLocalBasis DuneLocalBasis>
  template <typename Derived1, typename Derived2>
  void LocalBasis<DuneLocalBasis>::evaluateFunctionAndJacobian(const DomainType& local, Eigen::PlainObjectBase<Derived1>& N,
                                   Eigen::PlainObjectBase<Derived2>& dN) const {
    evaluateFunction(local, N);
    evaluateJacobian(local, dN);
  }

  template <Concepts::DuneLocalBasis DuneLocalBasis>
  template <typename IntegrationRule, typename... Ints>  requires std::conjunction_v<std::is_convertible<int, Ints>
                                ...>
      void LocalBasis<DuneLocalBasis>::bind(IntegrationRule&& p_rule, Impl::Derivatives<Ints...>&& ints) {
    rule             = p_rule;
    boundDerivatives = ints.args;
    Nbound           = std::make_optional(std::vector<Eigen::VectorX<RangeFieldType>>{});
    dNbound          = std::make_optional(std::vector<Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>>{});
    dNbound.value().resize(rule.value().size());
    Nbound.value().resize(rule.value().size());

    for (int i = 0; auto& gp : rule.value()) {
      if (boundDerivatives.value().contains(0)) evaluateFunction(gp.position(), Nbound.value()[i]);
      if (boundDerivatives.value().contains(1)) evaluateJacobian(gp.position(), dNbound.value()[i]);
      ++i;
    }
  }

}