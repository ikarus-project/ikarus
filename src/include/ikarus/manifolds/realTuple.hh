//
// Created by Alex on 19.05.2021.
//

#pragma once
#include <concepts>

#include <Eigen/Core>

namespace Ikarus::Manifold {
  /**
   * \brief FunctionReturnType of Euklidean space \f$\mathbb{R}^d\f$
   *
   * \tparam ct The type used for the scalar coordinate values, e.g. double, float
   * \tparam d Dimension of the embedding space of the manifold
   */

  template <typename ct, int d>
  class RealTuple {
  public:
    /** \brief Type used for coordinates */
    using ctype      = ct;
    using field_type = ct;

    /** \brief Size of how much values are needed to store the manifold */
    static constexpr int valueSize = d;

    /** \brief Size of how much values are needed to store the correction vector */
    static constexpr int correctionSize = d;

    /** \brief VectorType of the values of the manifold */
    using CoordinateType = Eigen::Matrix<ctype, valueSize, 1>;

    /** \brief VectorType of the values of the correction living in the tangentspace */
    using CorrectionType = Eigen::Matrix<ctype, correctionSize, 1>;

    RealTuple() = default;

    template <typename OtherType>
    struct Rebind {
      using type = RealTuple<OtherType, valueSize>;
    };

    /** \brief Copy-Constructor from the values in terms of coordinateType */
    explicit RealTuple(const CoordinateType &vec) noexcept : var{vec} {}

    /** \brief Move-Constructor from the values in terms of coordinateType */
    explicit RealTuple(CoordinateType &&vec) noexcept : var{std::move(vec)} {}

    /** \brief Get value of the manifold coordinates */
    CoordinateType getValue() const { return var; }

    /** \brief Set the coordinates of the manifold by const reference */
    void setValue(const CoordinateType &vec) { var = vec; }

    /** \brief Set the coordinates of the manifold by r_value reference */
    void setValue(CoordinateType &&vec) noexcept { var = std::move(vec); }

    /** \brief Update the manifold by an correction vector of size correctionSize */
    void update(const CorrectionType &correction) noexcept { var += correction; }

    /** \brief Access to data by const reference */
    const ctype &operator[](int i) const { return var[i]; }

    /** \brief Access to data by const reference */
    ctype &operator[](int i) { return var[i]; }

    auto &operator+=(const CorrectionType &correction) {
      this->update(correction);
      return *this;
    }

    /** \brief size */
    [[nodiscard]] size_t size() const { return var.size(); }
    auto begin() { return var.begin(); }
    auto end() { return var.end(); }

    auto begin() const { return var.begin(); }
    auto end() const { return var.end(); }

  private:
    CoordinateType var{CoordinateType::Zero()};
  };

  template <typename ctype2, int d2>
  std::ostream &operator<<(std::ostream &s, const RealTuple<ctype2, d2> &var2) {
    s << var2.getValue();
    return s;
  }

  template <typename ctype2, int d2, typename CorrectionType>
  [[nodiscard]] RealTuple<ctype2, d2> update(const RealTuple<ctype2, d2> &rt, const CorrectionType &correction) {
    return RealTuple<ctype2, d2>(rt.getValue() + correction);
  }

}  // namespace Ikarus::Manifold
