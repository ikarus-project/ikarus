//
// Created by Alex on 25.05.2021.
//

#include <ikarus/Variables/GenericVariable.h>

namespace Ikarus::Variable {

  int valueSize(const GenericVariable& vo) { return vo.variableImpl->do_valueSize(); }
  int correctionSize(const GenericVariable& vo) { return vo.variableImpl->do_correctionSize(); }
  void update(GenericVariable& vo, const Eigen::Ref<const Ikarus::DynVectord>& correction) {
    return vo.variableImpl->do_update(correction);
  }
  void setValue(GenericVariable& vo, const Eigen::Ref<const Ikarus::DynVectord>& value) {
    return vo.variableImpl->do_setValue(value);
  }
  Ikarus::DynVectord getValue(const GenericVariable& vo) { return vo.variableImpl->do_getValue(); }
  size_t getTag(const GenericVariable& var) { return var.variableImpl->do_getTag(); }
  bool operator==(const GenericVariable& var, const GenericVariable& other) {
    return var.variableImpl->do_equalComparison(other);
  }
  bool operator<(const GenericVariable& var, const GenericVariable& other) {
    return var.variableImpl->do_lessComparison(other);
  }

  std::ostream& operator<<(std::ostream& s, const GenericVariable& var) {
    s << var.variableImpl->do_getValue() << '\n'
      << " Tag: " << var.variableImpl->do_getTag() << '\n';
    return s;
  }
  size_t correctionSize(std::span<const GenericVariable> varSpan) {
    return std::accumulate(
        varSpan.begin(), varSpan.end(), 0,
        [](size_t cursize, const GenericVariable& var) { return cursize + correctionSize(var); });
  }
  void update(std::span<GenericVariable> varSpan, const Ikarus::DynVectord& correction) {
    assert(static_cast<Eigen::Index>(correctionSize(varSpan)) == correction.size());
    // update Variable
    Eigen::Index posHelper = 0;
    std::for_each(varSpan.begin(), varSpan.end(), [&](GenericVariable& var) {
      update(var, correction.segment(posHelper, correctionSize(var)));
      posHelper += correctionSize(var);
    });
  }

  size_t valueSize(std::span<const GenericVariable> varSpan) {
    return std::accumulate(
        varSpan.begin(), varSpan.end(), 0,
        [](size_t cursize, const GenericVariable& var) { return cursize + valueSize(var); });
  }

}  // namespace Ikarus::Variable