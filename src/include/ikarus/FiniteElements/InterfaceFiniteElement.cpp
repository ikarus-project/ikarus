//
// Created by Alex on 03.07.2021.
//

#include <map>
#include <utility>

#include <ikarus/FiniteElements/InterfaceFiniteElement.h>

namespace Ikarus::FiniteElements {

  void initialize(IFiniteElement& fe) { fe.feimpl->do_initialize(); }
  int dofSize(const IFiniteElement& fe) { return fe.feimpl->do_dofSize(); }
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const IFiniteElement& fe,
                                                                   const MatrixAffordances& matA,
                                                                   const VectorAffordances& vecA,
                                                                   IFiniteElement::VariableVectorType& vars,
                                                                   IFiniteElement::DataVectorType data) {
    return fe.feimpl->do_calculateLocalSystem(matA, vecA, vars, data);
  }
  Eigen::MatrixXd calculateMatrix(const IFiniteElement& fe, const MatrixAffordances& matA,
                                  IFiniteElement::VariableVectorType& vars, IFiniteElement::DataVectorType data) {
    return fe.feimpl->do_calculateMatrix(matA, vars, data);
  }

  Eigen::VectorXd calculateVector(const IFiniteElement& fe, const VectorAffordances& vecA,
                                  IFiniteElement::VariableVectorType& vars, IFiniteElement::DataVectorType data) {
    return fe.feimpl->do_calculateVector(vecA, vars, data);
  }

  double calculateScalar(const IFiniteElement& fe, const ScalarAffordances& scalA,
                         IFiniteElement::VariableVectorType& vars, IFiniteElement::DataVectorType data) {
    return fe.feimpl->do_calculateScalar(scalA, vars, data);
  }

  IFiniteElement::DofPairVectorType getEntityVariableTuple(const IFiniteElement& fe) {
    return fe.feimpl->do_getEntityVariableTuple();
  }
  unsigned int subEntities(const IFiniteElement& fe, unsigned int codim) { return fe.feimpl->do_subEntities(codim); }
  size_t subIndex(const IFiniteElement& fe, int i, unsigned int codim) { return fe.feimpl->do_subIndex(i, codim); }
  unsigned int dimension(const IFiniteElement& fe) { return fe.feimpl->do_dimension(); }
}  // namespace Ikarus::FiniteElements
