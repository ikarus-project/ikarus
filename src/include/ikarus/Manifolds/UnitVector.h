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
  template <typename ct, int d>
  class UnitVector {
  public:
    /** \brief Type used for coordinates */
    using ctype      = ct;
    using field_type = ct;

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

    /** \brief Access to data by const reference */
    const ctype &operator[](int i) const { return var[i]; }

    /** \brief Access to data by const reference */
    ctype &operator[](int i) { return var[i]; }

    /** \brief size */
    size_t size() const { return var.size(); }

    auto begin() { return var.begin(); }
    auto end() { return var.end(); }

    auto begin() const { return var.begin(); }
    auto end() const { return var.end(); }

    /** \brief Update the manifold by an correction vector of size correctionSize
     * For the unit vector in R^3 the correction are of size 2
     * Therefore, we need an basis for the tangent space.
     * This means we have two three dimensional vectors spanning this space.
     * This is done using the function orthonormalFrame which returns a 3x2 Matrix */
    void update(const CorrectionType &correction) {
      var += orthonormalFrame() * correction;
      var.normalize();  // projection-based retraction
    }

    static Eigen::Matrix<ctype, valueSize, valueSize> derivativeOfProjectionWRTposition(
        const Eigen::Vector<ctype, valueSize> &p) {
      const ctype normSquared = p.squaredNorm();
      const ctype norm        = std::sqrt(normSquared);

      Eigen::Matrix<ctype, valueSize, valueSize> result
          = (Eigen::Matrix<ctype, valueSize, valueSize>::Identity() - (p * p.transpose()) / normSquared) / norm;

      return result;
    }

    static Eigen::Matrix<ctype, valueSize, valueSize> secondDerivativeOfProjectionWRTposition(
        const Eigen::Vector<ctype, valueSize> &p, const Eigen::Ref<const Eigen::Vector<ctype, valueSize>> &along) {
      const ctype normSquared                  = p.squaredNorm();
      const ctype norm                         = std::sqrt(normSquared);
      const Eigen::Vector<ctype, valueSize> pN = p / norm;

      Eigen::Matrix<ctype, valueSize, valueSize> Q_along
          = 1 / normSquared
            * (pN.dot(along) * (3 * pN * pN.transpose() - Eigen::Matrix<ctype, valueSize, valueSize>::Identity())
               - along * pN.transpose() - pN * along.transpose());

      return Q_along;
    }

    static Eigen::Matrix<ctype, valueSize, valueSize> thirdDerivativeOfProjectionWRTposition(
        const Eigen::Vector<ctype, valueSize> &p,
        const Eigen::Ref<const Eigen::Vector<ctype, valueSize>> &along1,
        const Eigen::Ref<const Eigen::Vector<ctype, valueSize>> &along2) {
      using FieldMat                           = Eigen::Matrix<ctype, valueSize, valueSize>;
      const ctype normSquared                  = p.squaredNorm();
      const ctype norm                         = std::sqrt(normSquared);
      const Eigen::Vector<ctype, valueSize> pN = p / norm;
      const ctype tscala1                      = pN.dot(along1);
      const ctype tscalwd1                     = pN.dot(along2);
      const ctype a1scalwd1                    = along1.dot(along2);
      const ctype normwcubinv                  = 1 / (normSquared * norm);
      const FieldMat a1dyadt                   = along1 * pN.transpose();
      const FieldMat wd1dyadt                  = along2 * pN.transpose();
      const FieldMat tDyadict                  = pN * pN.transpose();
      const FieldMat Id3minus5tdyadt           = FieldMat::Identity() - 5.0 * tDyadict;
      FieldMat Chi_along
          = normwcubinv
            * (3.0 * tscalwd1 * (a1dyadt + 0.5 * tscala1 * Id3minus5tdyadt)
               + 3.0 * (0.5 * a1scalwd1 * tDyadict + tscala1 * wd1dyadt) - along1 * along2.transpose()
               - a1scalwd1 * 0.5 * FieldMat::Identity());
      Chi_along = (Chi_along + Chi_along.transpose()).eval();
      return Chi_along;
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

    auto &operator+=(const CorrectionType &correction) {
      this->update(correction);
      return *this;
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

  template <typename ctype2, int d2>
  [[nodiscard]] UnitVector<ctype2, d2> operator+(const UnitVector<ctype2, d2> &rt,
                                                 const typename UnitVector<ctype2, d2>::CorrectionType &correction) {
    return UnitVector<ctype2, d2>(rt.getValue() + rt.orthonormalFrame() * correction);
  }

}  // namespace Ikarus::Manifold
