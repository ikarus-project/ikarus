//
// Created by Alex on 03.07.2021.
//

#include <map>
#include <utility>

#include <ikarus/FiniteElements/InterfaceFiniteElement.h>

namespace Ikarus::FiniteElements {

  int dofSize(const IFiniteElement& fe) { return fe.feimpl->do_dofSize(); }
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const IFiniteElement& fe,
                                                                   const IFiniteElement::FEParameterType & par) {
    return fe.feimpl->do_calculateLocalSystem(par);
  }
  Eigen::MatrixXd calculateMatrix(const IFiniteElement& fe, const IFiniteElement::FEParameterType & par) {
    return fe.feimpl->do_calculateMatrix(par);
  }

  Eigen::VectorXd calculateVector(const IFiniteElement& fe, const IFiniteElement::FEParameterType & par) {
    return fe.feimpl->do_calculateVector(par);
  }

  double calculateScalar(const IFiniteElement& fe, const IFiniteElement::FEParameterType & par) {
    return fe.feimpl->do_calculateScalar(par);
  }

  IFiniteElement::DofPairVectorType getEntityVariableTuple(const IFiniteElement& fe) {
    return fe.feimpl->do_getEntityVariableTuple();
  }
  unsigned int subEntities(const IFiniteElement& fe, unsigned int codim) { return fe.feimpl->do_subEntities(codim); }
  size_t subIndex(const IFiniteElement& fe, int i, unsigned int codim) { return fe.feimpl->do_subIndex(i, codim); }
  unsigned int dimension(const IFiniteElement& fe) { return fe.feimpl->do_dimension(); }
}  // namespace Ikarus::FiniteElements
