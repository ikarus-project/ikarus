//
// Created by Alex on 03.07.2021.
//

#include <map>
#include <utility>

#include "InterfaceFiniteElement.h"

template <typename LocalView, typename SolutionType>
Eigen::MatrixXd IFiniteElement<LocalView,SolutionType>::calculateMatrix(
    const IFiniteElement<LocalView,SolutionType>::FERequirementType& par) const {
  return feimpl->do_calculateMatrix(par);
}
template <typename LocalView, typename SolutionType>
Eigen::VectorXd IFiniteElement<LocalView,SolutionType>::calculateVector(
    const IFiniteElement<LocalView,SolutionType>::FERequirementType& par) const {
  return feimpl->do_calculateVector(par);
}

template <typename LocalView, typename SolutionType>
double IFiniteElement<LocalView,SolutionType>::calculateScalar(const IFiniteElement<LocalView,SolutionType>::FERequirementType& par) const {
  return feimpl->do_calculateScalar(par);
}


template <typename LocalView, typename SolutionType>
std::vector<typename IFiniteElement<LocalView,SolutionType>::GlobalIndex> IFiniteElement<LocalView,SolutionType>::globalIndices() const {
  return feimpl->do_globalIndices();
}

template <typename LocalView, typename SolutionType>
std::pair<Eigen::MatrixXd, Eigen::VectorXd> IFiniteElement<LocalView,SolutionType>::calculateLocalSystem( const IFiniteElement<LocalView,SolutionType>::FERequirementType &par) const {
  return feimpl->do_calculateLocalSystem();
}