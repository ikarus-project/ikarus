//
// Created by Alex on 25.05.2021.
//

#include <ikarus/Variables/InterfaceVariable.h>
#include <ikarus/Variables/VariableDefinitions.h>

namespace Ikarus::Variable {

  int valueSize(const IVariable& vo) { return vo.variableImpl->do_valueSize(); }
  int correctionSize(const IVariable& vo) { return vo.variableImpl->do_correctionSize(); }

  IVariable& operator-=(IVariable& vo, const IVariable::UpdateType& correction) {
    vo.variableImpl->do_assignAdd(-correction);
    return vo;
  }

  IVariable& operator+=(IVariable& vo, const IVariable::UpdateType& correction) {
    vo.variableImpl->do_assignAdd(correction);
    return vo;
  }

  IVariable& operator+=(IVariable* vo, const IVariable::UpdateType& correction) { return ((*vo) += correction); }
  IVariable& operator-=(IVariable* vo, const IVariable::UpdateType& correction) { return ((*vo) -= correction); }

  IVariable operator+(IVariable& vo, const IVariable::UpdateType& correction) {
    IVariable res{vo};
    res.variableImpl->do_assignAdd(correction);
    return res;
  }

  IVariable operator-(IVariable& vo, const IVariable::UpdateType& correction) {
    IVariable res{vo};
    res.variableImpl->do_assignAdd(-correction);
    return res;
  }

  IVariable operator+(IVariable* vo, const IVariable::UpdateType& correction) { return ((*vo) + correction); }
  IVariable operator-(IVariable* vo, const IVariable::UpdateType& correction) { return ((*vo) - correction); }

  void setValue(IVariable& vo, const IVariable::UpdateType& value) { return vo.variableImpl->do_setValue(value); }
  IVariable::CoordinateType getValue(const IVariable& vo) { return vo.variableImpl->do_getValue(); }
  IVariable::CoordinateType getValue(const IVariable* vo) { return getValue(*vo); }
  int getTag(const IVariable& var) { return var.variableImpl->do_getTag(); }
  bool operator==(const IVariable& var, const IVariable& other) { return var.variableImpl->do_equalComparison(other); }
  bool operator<(const IVariable& var, const IVariable& other) { return var.variableImpl->do_lessComparison(other); }

  std::ostream& operator<<(std::ostream& s, const IVariable& var) {
    s << var.variableImpl->do_getValue().transpose() << '\n' << " Tag: " << getName(var) << '\n';
    return s;
  }

  std::ostream& operator<<(std::ostream& s, const IVariable* var) {
    s << (*var);
    return s;
  }

  size_t correctionSize(std::span<const IVariable> varSpan) {
    return std::accumulate(varSpan.begin(), varSpan.end(), size_t{0},
                           [](size_t cursize, const IVariable& var) { return cursize + correctionSize(var); });
  }
  void update(std::span<IVariable> varSpan, const Eigen::VectorXd& correction) {
    assert(static_cast<Eigen::Index>(correctionSize(varSpan)) == correction.size());
    // update Variable
    Eigen::Index posHelper = 0;
    std::for_each(varSpan.begin(), varSpan.end(), [&](IVariable& var) {
      var += correction.segment(posHelper, correctionSize(var));
      posHelper += correctionSize(var);
    });
  }

  size_t valueSize(std::span<const IVariable> varSpan) {
    return std::accumulate(varSpan.begin(), varSpan.end(), size_t{0},
                           [](size_t cursize, const IVariable& var) { return cursize + valueSize(var); });
  }

  std::string getName(const IVariable& var) { return Ikarus::Variable::variableNames[getTag(var)]; }

  bool isType(const IVariable& vo, Ikarus::Variable::VariableTags tag) { return getTag(vo) == static_cast<int>(tag); }
  bool isType(IVariable* vo, Ikarus::Variable::VariableTags tag) { return isType(*vo, tag); }

  double& IVariable::operator[](int i) { return this->variableImpl->operator[](i); }

}  // namespace Ikarus::Variable