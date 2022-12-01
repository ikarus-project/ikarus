// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

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


  /*
   * This function returns the second derivatives of the ansatz functions.
   * The assumed order is in Voigt notation, e.g. for 3d ansatzfunctions N_xx,N_yy,N_zz,N_yz,N_xz, N_xy
   */
  template <Concepts::DuneLocalBasis DuneLocalBasis>
  template <typename Derived>
  void LocalBasis<DuneLocalBasis>::evaluateSecondDerivatives(const DomainType& local, Eigen::PlainObjectBase<Derived>& ddN) const {
    std::array<unsigned int, gridDim> order;
    std::ranges::fill(order, 0);
    ddN.setZero(dNdune.size(), Eigen::NoChange);

    for (int i = 0; i < gridDim; ++i) { //Diagonal terms
      order[i] = 2;
      duneLocalBasis->partial(order,local, ddNdune);
      for (size_t j = 0; j < ddNdune.size(); ++j)
        ddN(j,i)=ddNdune[j][0];

      order[i] = 0;
    }

    std::ranges::fill(order, 1);
    for (int i = 0; i < gridDim*(gridDim-1)/2; ++i) { //off-diagonal terms
      if constexpr (gridDim>2)
        order[i] = 0;
      duneLocalBasis->partial(order,local, ddNdune);
      for (size_t j = 0; j < ddNdune.size(); ++j)
        ddN(j,i+gridDim)=ddNdune[j][0];
      order[i] = 1;
    }

  }

  template <Concepts::DuneLocalBasis DuneLocalBasis>
  template <typename Derived1, typename Derived2>
  void LocalBasis<DuneLocalBasis>::evaluateFunctionAndJacobian(const DomainType& local, Eigen::PlainObjectBase<Derived1>& N,
                                   Eigen::PlainObjectBase<Derived2>& dN) const {
    evaluateFunction(local, N);
    evaluateJacobian(local, dN);
  }

  template <Concepts::DuneLocalBasis DuneLocalBasis>
  template <typename... Ints>  requires std::conjunction_v<std::is_convertible<int, Ints>
                                ...>
      void LocalBasis<DuneLocalBasis>::bind(const Dune::QuadratureRule<DomainFieldType, gridDim>& p_rule, Impl::Derivatives<Ints...>&& ints) {
    rule             = p_rule;
    boundDerivatives = ints.args;
    Nbound           = std::make_optional<typename decltype(Nbound)::value_type> ();
    dNbound           = std::make_optional<typename decltype(dNbound)::value_type> ();
    ddNbound           = std::make_optional<typename decltype(ddNbound)::value_type> ();
    dNbound.value().resize(rule.value().size());
    ddNbound.value().resize(rule.value().size());
    Nbound.value().resize(rule.value().size());

    for (int i = 0; auto& gp : rule.value()) {
      if (boundDerivatives.value().contains(0)) evaluateFunction(gp.position(), Nbound.value()[i]);
      if (boundDerivatives.value().contains(1)) evaluateJacobian(gp.position(), dNbound.value()[i]);
      if (boundDerivatives.value().contains(2)) evaluateSecondDerivatives(gp.position(), ddNbound.value()[i]);
      ++i;
    }
  }

}
