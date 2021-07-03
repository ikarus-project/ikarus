//
// Created by Alex on 19.05.2021.
//

#pragma once
#include <concepts>

#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::Manifold {
  /**
   * \brief Manifold of unit vectors \f$\mathcal{S}^{d-1}\f$ embedded into space \f$\mathbb{R}^d\f$
   *
   * \tparam ct The type used for the scalar coordinate values, e.g. double,float
   * \tparam d Dimension of the embedding space of the manifold
   */
  template <std::floating_point ct, int d>
  class UnitVector {
  public:
    /** \brief Type used for coordinates */
    using ctype = ct;

    /** \brief Size of how much values are needed to store the manifold */
    static constexpr int valueSize = d;

    /** \brief Size of how much values are needed to store the correction vector */
    static constexpr int correctionSize = d - 1;

    /** \brief VectorType of the values of the manifold */
    using CoordinateType = typename Eigen::Vector<ctype, valueSize>;

    /** \brief VectorType of the values of the correction living in the tangentspace */
    using CorrectionType = typename Eigen::Vector<ctype, correctionSize>;

    UnitVector() = default;

    ~UnitVector()                  = default;                 // destructor
    UnitVector(const UnitVector &) = default;                 // copy constructor
    UnitVector &operator=(const UnitVector &) = default;      // copy assignment
    UnitVector(UnitVector &&) noexcept        = default;      // move constructor
    UnitVector &operator=(UnitVector &&) noexcept = default;  // move assignment

    /** \brief Copy-Constructor from the values in terms of coordinateType */
    explicit UnitVector(const CoordinateType &vec) noexcept : var{vec.normalized()} {}

    /** \brief Move-Constructor from the values in terms of coordinateType */
    explicit UnitVector(CoordinateType &&vec) noexcept : var{vec.normalized()} {}

    /** \brief Get the coordinates of the manifold by value */
    CoordinateType getValue() const { return var; }

    /** \brief Set the coordinates of the manifold by const reference */
    void setValue(const CoordinateType &vec) { var = vec.normalized(); }

    /** \brief Set the coordinates of the manifold by r_value reference */
    void setValue(CoordinateType &&vec) { var = std::move(vec.normalized()); }

    /** \brief Update the manifold by an correction vector of size correctionSize
     * For the unit vector in R^3 the correction are of size 2
     * Therefore, we need an basis for the tangent space.
     * This means we have two three dimensional vectors spanning this space.
     * This is done using the function orthonormalFrame which returns a 3x2 Matrix */
    void update(const CorrectionType &correction) {
      var += orthonormalFrame() * correction;
      var.normalize();  // projection-based retraction
    }

    /** \brief Compute an orthonormal basis of the tangent space of S^n.
     * Taken from Oliver Sander dune-gfe */
    Eigen::Matrix<ctype, valueSize, correctionSize> orthonormalFrame() const {
      Eigen::Matrix<ctype, valueSize, correctionSize> result;

      // Coordinates of the stereographic projection
      Eigen::Matrix<ctype, correctionSize, 1> X;

      if (var[valueSize - 1] <= 0) {
        // Stereographic projection from the north pole onto R^{N-1}
        for (size_t i = 0; i < valueSize - 1; i++)
          X[i] = var[i] / (1 - var[valueSize - 1]);

      } else {
        // Stereographic projection from the south pole onto R^{N-1}
        for (size_t i = 0; i < valueSize - 1; i++)
          X[i] = var[i] / (1 + var[valueSize - 1]);
      }

      ctype RSquared = X.squaredNorm();

      for (size_t i = 0; i < valueSize - 1; i++)
        for (size_t j = 0; j < valueSize - 1; j++)
          // Note: the matrix is the transpose of the one in the paper
          result(j, i) = 2 * (i == j) * (1 + RSquared) - 4 * X[i] * X[j];

      for (size_t j = 0; j < valueSize - 1; j++)
        result(valueSize - 1, j) = 4 * X[j];

      // Upper hemisphere: adapt formulas so it is the stereographic projection from the south pole
      if (var[valueSize - 1] > 0)
        for (size_t j = 0; j < valueSize - 1; j++)
          result(valueSize - 1, j) *= -1;

      // normalize the rows to make the orthogonal basis orthonormal
      for (size_t i = 0; i < valueSize - 1; i++)
        result.col(i).normalize();

      return result;
    }

  private:
    CoordinateType var{CoordinateType::UnitX()};
  };

  template <typename ctype2, int d2>
  std::ostream &operator<<(std::ostream &s, const UnitVector<ctype2, d2> &var2) {
    s << var2.getValue();
    return s;
  }

  template <typename ctype2, int d2>
  [[nodiscard]] UnitVector<ctype2, d2> update(const UnitVector<ctype2, d2> &rt,
                                              const typename UnitVector<ctype2, d2>::CorrectionType &correction) {
    return UnitVector<ctype2, d2>(rt.getValue() + rt.orthonormalFrame() * correction);
  }

}  // namespace Ikarus::Manifold
