//
// Created by Alex on 19.05.2021.
//

#pragma once

#include <variant>

#include "ikarus/Manifolds/ManifoldInterface.h"

namespace Ikarus {

  template <Concepts::Manifold Mani, int tag>
  class DefaultVariable {
  public:
    /** \brief Type of the impl manifold, usually Realtuple */
    using ManifoldType = Mani;

    /** \brief Size of how much values are needed to store the variable */
    static constexpr int valueSize = ManifoldType::valueSize;

    /** \brief Size of how much values are needed to store the correction vector */
    static constexpr int correctionSize = ManifoldType::correctionSize;

    /** \brief Type used for coordinates */
    using ctype = typename ManifoldType::ctype;

    /** \brief VectorType of the values of the variable */
    using CoordinateType = typename ManifoldType::CoordinateType;

    /** \brief VectorType of the values of the correction */
    using CorrectionType = typename ManifoldType::CorrectionType;

    /** \brief Value of the tag to distinguish variables  */
    static constexpr auto tagvalue = tag;

    DefaultVariable()                       = default;
    ~DefaultVariable()                      = default;                 // destructor
    DefaultVariable(const DefaultVariable&) = default;                 // copy constructor
    DefaultVariable& operator=(const DefaultVariable&) = default;      // copy assignment
    DefaultVariable(DefaultVariable&&) noexcept        = default;      // move constructor
    DefaultVariable& operator=(DefaultVariable&&) noexcept = default;  // move assignment

    /** \brief Copy-Constructor from the values in terms of coordinateType */
    explicit DefaultVariable(const CoordinateType& vec) noexcept : var{vec} {}

    /** \brief Move-Constructor from the values in terms of coordinateType */
    explicit DefaultVariable(CoordinateType&& vec) noexcept : var{std::move(vec)} {}

    /** \brief Get the value of the variable coordinates */
    CoordinateType getValue() const { return var.getValue(); }

    /** \brief Set the coordinates of the variable by const reference */
    void setValue(const CoordinateType& vec) { var.setValue(vec); }

    /** \brief Set the coordinates of the variable by r_value reference */
    void setValue(CoordinateType&& vec) noexcept { var.setValue(std::move(vec)); }

    /** \brief Returns the tag of variable */
    [[nodiscard]] static constexpr auto getTag() { return tagvalue; }

    void update(const CorrectionType& correction) { var.update(correction); }

    /** \brief Access to data by const reference */
    const ctype& operator[](int i) const { return var[i]; }

    /** \brief Access to data by const reference */
    ctype& operator[](int i) { return var[i]; }

  private:
    ManifoldType var{};
  };

  template <typename Mani, int tag>
  [[nodiscard]] DefaultVariable<Mani, tag> update(
      const DefaultVariable<Mani, tag>& rt, const typename DefaultVariable<Mani, tag>::CorrectionType& correction) {
    auto res{rt};
    res.update(correction);
    return res;
  }

  template <typename Mani, int tag, typename Mani2, int tag2>
  constexpr inline bool operator==(const DefaultVariable<Mani, tag>& lhs, const DefaultVariable<Mani2, tag2>& rhs) {
    return lhs.getTag() == rhs.getTag();
  }

  template <typename Mani, int tag, typename Mani2, int tag2>
  constexpr inline bool operator!=(const DefaultVariable<Mani, tag>& lhs, const DefaultVariable<Mani2, tag2>& rhs) {
    return !(lhs == rhs);
  }

  template <typename Mani, int tag>
  std::ostream& operator<<(std::ostream& s, const DefaultVariable<Mani, tag>& var) {
    s << var.getValue().transpose() << "\n";
    return s;
  }

}  // namespace Ikarus
