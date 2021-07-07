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
                                                                   IFiniteElement::VariableVectorType& vars,
                                                                   const MatrixAffordances& matA,
                                                                   const VectorAffordances& vecA) {
    return fe.feimpl->do_calculateLocalSystem(vars, matA, vecA);
  }
  Eigen::MatrixXd calculateMatrix(const IFiniteElement& fe, IFiniteElement::VariableVectorType& vars,
                                  const MatrixAffordances& matA) {
    return fe.feimpl->do_calculateMatrix(vars, matA);
  }

  Eigen::MatrixXd calculateMatrix(const IFiniteElement* fe, IFiniteElement::VariableVectorType& vars,
                                  const MatrixAffordances& matA) {
    return calculateMatrix((*fe), vars, matA);
  }

  Eigen::VectorXd calculateVector(const IFiniteElement& fe, IFiniteElement::VariableVectorType& vars,
                                  const VectorAffordances& vecA) {
    return fe.feimpl->do_calculateVector(vars, vecA);
  }
  Eigen::VectorXd calculateVector(const IFiniteElement* fe, IFiniteElement::VariableVectorType& vars,
                                  const VectorAffordances& vecA) {
    return calculateVector((*fe), vars, vecA);
  }

  double calculateScalar(const IFiniteElement& fe, IFiniteElement::VariableVectorType& vars,
                         const ScalarAffordances& scalA) {
    return fe.feimpl->do_calculateScalar(vars, scalA);
  }

  double calculateScalar(const IFiniteElement* fe, IFiniteElement::VariableVectorType& vars,
                         const ScalarAffordances& scalA) {
    return calculateScalar((*fe),vars,scalA);
  }
  IFiniteElement::DofPairVectorType getEntityVariablePairs(const IFiniteElement& fe) {
    return fe.feimpl->do_getEntityVariablePairs();
  }
  size_t getEntityID(const IFiniteElement& fe) { return fe.feimpl->do_getEntityID(); }
}  // namespace Ikarus::FiniteElements


