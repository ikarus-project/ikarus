//
// Created by lex on 04/02/2022.
//

#pragma once
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <Eigen/Core>

namespace Ikarus {
  template <typename DuneLocalBasis>
  class LocalBasis {
  public:
    LocalBasis(const DuneLocalBasis& p_basis) : duneLocalBasis{&p_basis} {}

    static constexpr int gridDim = DuneLocalBasis::Traits::dimDomain;
    using DomainType             = typename DuneLocalBasis::Traits::DomainType;
    using RangeType              = typename DuneLocalBasis::Traits::RangeType;
    using DomainFieldType        = typename DuneLocalBasis::Traits::DomainFieldType;
    using RangeFieldType         = typename DuneLocalBasis::Traits::RangeFieldType;
    using JacobianType           = typename DuneLocalBasis::Traits::JacobianType;

    void evaluateFunction(const DomainType& local, Eigen::VectorX<RangeFieldType>& N) {
      duneLocalBasis->evaluateFunction(local, Ndune);
      N.setZero(Ndune.size());
      for (size_t i = 0; i < Ndune.size(); ++i)
        N[i] = Ndune[i][0];
    }

    void evaluateJacobian(const DomainType& local, Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& dN) {
      duneLocalBasis->evaluateJacobian(local, dNdune);
      dN.setZero(dNdune.size(), Eigen::NoChange);
      for (auto i = 0U; i < dNdune.size(); ++i)
        for (int j = 0; j < gridDim; ++j)
          dN(i, j) = dNdune[i][0][j];
    }

    void evaluateFunctionAndJacobian(const DomainType& local, Eigen::VectorX<RangeFieldType>& N,
                                     Eigen::Matrix<RangeFieldType, Eigen::Dynamic, gridDim>& dN) {
      evaluateFunction(local, N);
      evaluateJacobian(local, dN);
    }

  private:
    std::vector<JacobianType> dNdune;
    std::vector<RangeType> Ndune;
    DuneLocalBasis const* duneLocalBasis;
  };

}  // namespace Ikarus