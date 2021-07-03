//
// Created by Alex on 03.07.2021.
//

#include <ikarus/FiniteElements/InterfaceFiniteElement.h>

namespace Ikarus::FiniteElements {
  void initialize(IFiniteElement& fe) { fe.feimpl->do_initialize(); }
  int dofSize(const IFiniteElement& fe) { return fe.feimpl->do_dofSize(); }
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculateLocalSystem(const IFiniteElement& fe,
                                                                   const ElementMatrixAffordances& matA,
                                                                   const ElementVectorAffordances& vecA) {
    return fe.feimpl->do_calculateLocalSystem(matA, vecA);
  }
  Eigen::MatrixXd calculateMatrix(const IFiniteElement& fe, const ElementMatrixAffordances& matA) {
    return fe.feimpl->do_calculateMatrix(matA);
  }
  Eigen::VectorXd calculateVector(const IFiniteElement& fe, const ElementVectorAffordances& vecA) {
    return fe.feimpl->do_calculateVector(vecA);
  }
  double calculateScalar(const IFiniteElement& fe, const ElementScalarAffordances& scalA) {
    return fe.feimpl->do_calculateScalar(scalA);
  }
  IFiniteElement::DofVectorType getEntityVariablePairs(const IFiniteElement& fe) {
    return fe.feimpl->do_getEntityVariablePairs();
  }
}  // namespace Ikarus::FiniteElements