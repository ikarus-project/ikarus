//
// Created by Alex on 19.05.2021.
//

#pragma once
#include <concepts>

#include <ikarus/utils/eigenDuneTransformations.hh>

namespace Ikarus {

  template<typename Scalar>
  Scalar arcCosSquared(const Scalar& x) {
    using std::acos;
    const Scalar eps = 1e-2;
    if (x > 1-eps) {  // acos is not differentiable, use the series expansion instead,
      // we need here lots of terms to be sure that the numerical derivatives are also within maschine precission
      return 11665028.0/4729725.0
             -141088.0/45045.0*x
             +   413.0/429.0*x*x
             -  5344.0/12285.0*Dune::power(x,3)
             +    245.0/1287.0*Dune::power(x,4)
             -  1632.0/25025.0*Dune::power(x,5)
             +     56.0/3861.0*Dune::power(x,6)
             -    32.0/21021.0*Dune::power(x,7);
    } else {
      return Dune::power(acos(x),2);
    }
  }

  template<typename Scalar>
  Scalar derivativeOfArcCosSquared(const Scalar& x) {
    using std::acos;
    using std::sqrt;
    const Scalar eps = 1e-2;
    if (x > 1-eps) {  // regular expression is unstable, use the series expansion instead
      // we need here lots of terms to be sure that the numerical derivatives are also within maschine precission
      //return -2 + 2*(x-1)/3 - 4/15*(x-1)*(x-1);
      return -47104.0/15015.0
             +12614.0/6435.0*x
             -63488.0/45045.0*x*x
             + 1204.0/1287.0*Dune::power(x,3)
             - 2048.0/4095.0*Dune::power(x,4)
             +   112.0/585.0*Dune::power(x,5)
             -2048.0/45045.0*Dune::power(x,6)
             +   32.0/6435.0*Dune::power(x,7);

    } else if (x < -1+eps) {  // The function is not differentiable
      DUNE_THROW(Dune::Exception, "arccos^2 is not differentiable at x==-1!"<< x);
    } else
      return -2*acos(x) / sqrt(1-x*x);
  }


  /** \brief Compute the second derivative of arccos^2 without getting unstable for x close to 1 */
  template<typename Scalar>
   Scalar secondDerivativeOfArcCosSquared(const Scalar& x) {
    using std::acos;
    using std::pow;
    const Scalar eps = 1e-2;
    if (x > 1-eps) {  // regular expression is unstable, use the series expansion instead
      // we need here lots of terms to be sure that the numerical derivatives are also within maschine precission
      //return 2.0/3 - 8*(x-1)/15;
      return 1350030.0/676039.0+5632.0/2028117.0*Dune::power(x,10)
             -1039056896.0/334639305.0*x
             +150876.0/39767.0*x*x
             -445186048.0/111546435.0*Dune::power(x,3)
             +       343728.0/96577.0*Dune::power(x,4)
             -  57769984.0/22309287.0*Dune::power(x,5)
             +      710688.0/482885.0*Dune::power(x,6)
             -  41615360.0/66927861.0*Dune::power(x,7)
             +     616704.0/3380195.0*Dune::power(x,8)
             -     245760.0/7436429.0*Dune::power(x,9);
    } else if (x < -1+eps) {  // The function is not differentiable
      DUNE_THROW(Dune::Exception, "arccos^2 is not differentiable at x==-1!");
    } else
      return 2/(1-x*x) - 2*x*acos(x) / pow(1-x*x,1.5);
  }




  /**
   * \brief FunctionReturnType of unit vectors \f$\mathcal{S}^{d-1}\f$ embedded into space \f$\mathbb{R}^d\f$
   *
   * \tparam ct The type used for the scalar coordinate values, e.g. double,float
   * \tparam d Dimension of the embedding space of the manifold
   */
  template <typename ct, int d>
  class UnitVector {
  public:
    /** \brief Type used for coordinates */
    using ctype      = ct;
    using Scalar      = ct;
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

    /** \brief Copy-Constructor from the values in terms of coordinateType */
    explicit UnitVector(const CoordinateType &vec) noexcept : var{vec.normalized()} {}

    /** \brief Move-Constructor from the values in terms of coordinateType */
    explicit UnitVector(CoordinateType &&vec) noexcept : var{vec.normalized()} {}

    const CoordinateType &getValue() const { return var; }

    void setValue(const CoordinateType &vec) { var = vec.normalized(); }

    /** \brief Set the coordinates of the manifold by r_value reference */
    void setValue(CoordinateType &&vec) { var = std::move(vec.normalized()); }

    /** \brief Distance between two elments */
    template<typename Scalar>
     auto squaredDistance(const UnitVector<Scalar,d> &other) const {
      const auto dotprod = var.dot(other.var);
      using std::min;
      dotprod = min(dotprod,1.0);
  if(acos(dotprod)> std::numbers::pi/2)
    std::cout<<"LARGE angle detected: "<<acos(dotprod)<<std::endl;
       return arcCosSquared(dotprod);
     }

   private:
     template<typename Scalar>
     CoordinateType euclideanDerivativeOfDistanceSquaredWRTSecondArgument( const UnitVector<Scalar,d> & other) const {
       const auto dotprod = var.dot(other.var);

       return var*derivativeOfArcCosSquared(dotprod);
     }
   public:

    template<typename Scalar>
    CorrectionType derivativeOfDistanceSquaredWRTSecondArgument( const UnitVector<Scalar,d> & other) const {
      return other.orthonormalFrame().transpose()*euclideanDerivativeOfDistanceSquaredWRTSecondArgument(other);
    }

     template<typename Scalar>
     Eigen::Matrix<Scalar,d-1,d-1> secondDerivativeOfDistanceSquaredWRTSecondArgument( const UnitVector<Scalar,d> & other) const {

       const auto dotprod = var.dot(other.var);
       const auto BLA = other.orthonormalFrame();
       const auto gradEuk = euclideanDerivativeOfDistanceSquaredWRTSecondArgument(other);
       return (BLA.transpose()*secondDerivativeOfArcCosSquared(dotprod)*var*var.transpose()*BLA-other.getValue().dot(gradEuk)*Eigen::Matrix<double,correctionSize,correctionSize>::Identity()).eval();
     }



    /** \brief Access to data by const reference */
    const ctype &operator[](int i) const { return var[i]; }

    /** \brief Access to data by const reference */
    ctype &operator[](int i) { return var[i]; }

    /** \brief size */
    [[nodiscard]] size_t size() const { return var.size(); }

    auto begin() { return var.begin(); }
    auto end() { return var.end(); }

    auto begin() const { return var.begin(); }
    auto end() const { return var.end(); }

    template <typename ctOther, int dOther>
    requires std::convertible_to<ctOther, ctype>
    friend class UnitVector;

    /** \brief Copy assignement if the other type has different underlying type*/
    template <typename ctype_>
    requires std::convertible_to<ctype_, ctype>
    UnitVector<ctype, d>
    &operator=(const UnitVector<ctype_, d> &other) {
      var = other.var;
      return *this;
    }

    template <typename OtherType>
    struct Rebind {
      using other = UnitVector<OtherType, valueSize>;
    };

    CoordinateType projectOntoTangentSpace(const CoordinateType& vec)
    {
      return derivativeOfProjectionWRTposition(var)*vec;
    }



    /** \brief Update the manifold by an correction vector of size correctionSize
     * For the unit vector in R^3 the correction are of size 2
     * Therefore, we need an basis for the tangent space.
     * This means we have two three dimensional vectors spanning this space.
     * This is done using the function orthonormalFrame which returns a 3x2 Matrix */
    void update(const CorrectionType &correction) {
      var += orthonormalFrame() * correction;
      var.normalize();  // projection-based retraction
    }

    void addInEmbedding(const CoordinateType &correction) {
      var+=correction;
    }

    static Eigen::Matrix<ctype, valueSize, valueSize> derivativeOfProjectionWRTposition(
        const Eigen::Vector<ctype, valueSize> &p) {
      const ctype norm                         = p.norm();
      const Eigen::Vector<ctype, valueSize> pN = p / norm;

      Eigen::Matrix<ctype, valueSize, valueSize> result
          = (Eigen::Matrix<ctype, valueSize, valueSize>::Identity() - pN * pN.transpose()) / norm;

      return result;
    }

    template <typename Derived>
    static Eigen::Matrix<ctype, valueSize, valueSize> secondDerivativeOfProjectionWRTposition(
        const Eigen::Vector<ctype, valueSize> &p, const Eigen::MatrixBase<Derived> &along) {
      const ctype normSquared = p.squaredNorm();
      using std::sqrt;
      const ctype norm                         = sqrt(normSquared);
      const Eigen::Vector<ctype, valueSize> pN = p / norm;

      Eigen::Matrix<ctype, valueSize, valueSize> Q_along
          = 1 / normSquared
            * (pN.dot(along) * (3 * pN * pN.transpose() - Eigen::Matrix<ctype, valueSize, valueSize>::Identity())
               - along * pN.transpose() - pN * along.transpose());

      return Q_along;
    }

    static Eigen::Matrix<ctype, valueSize, valueSize> thirdDerivativeOfProjectionWRTposition(
        const Eigen::Vector<ctype, valueSize> &p, const Eigen::Ref<const Eigen::Vector<ctype, valueSize>> &along1,
        const Eigen::Ref<const Eigen::Vector<ctype, valueSize>> &along2) {
      using FieldMat          = Eigen::Matrix<ctype, valueSize, valueSize>;
      const ctype normSquared = p.squaredNorm();
      using std::sqrt;
      const ctype norm                         = sqrt(normSquared);
      const Eigen::Vector<ctype, valueSize> pN = p / norm;
      const ctype tscala1                      = pN.dot(along1);
      const ctype tscalwd1                     = pN.dot(along2);
      const ctype a1scalwd1                    = along1.dot(along2);
      const ctype normwcubinv                  = 1 / (normSquared * norm);
      const FieldMat a1dyadt                   = along1 * pN.transpose();
      const FieldMat wd1dyadt                  = along2 * pN.transpose();
      const FieldMat tDyadict                  = pN * pN.transpose();
      const FieldMat Id3minus5tdyadt           = FieldMat::Identity() - 5.0 * tDyadict;
      FieldMat Chi_along                       = normwcubinv
                           * (3.0 * tscalwd1 * (a1dyadt + 0.5 * tscala1 * Id3minus5tdyadt)
                              + 3.0 * (0.5 * a1scalwd1 * tDyadict + tscala1 * wd1dyadt) - along1 * along2.transpose()
                              - a1scalwd1 * 0.5 * FieldMat::Identity());
      Chi_along = (Chi_along + Chi_along.transpose()).eval();
      return Chi_along;
    }

    /** \brief Compute an orthonormal basis of the tangent space of S^n.
     * Taken from Oliver Sander's dune-gfe */
    Eigen::Matrix<ctype, valueSize, correctionSize> orthonormalFrame() const {
      using ResultType = Eigen::Matrix<ctype, valueSize, correctionSize>;
      ResultType result;

      // Coordinates of the stereographic projection
      Eigen::Matrix<ctype, correctionSize, 1> X;

      if (var[correctionSize] <= 0)
        // Stereographic projection from the north pole onto R^{N-1}
        X = var.template head<correctionSize>() / (1 - var[correctionSize]);
      else
        // Stereographic projection from the south pole onto R^{N-1}
        X = var.template head<correctionSize>() / (1 + var[correctionSize]);

      result.template topLeftCorner<correctionSize, correctionSize>()
          = (2 * (1 + X.squaredNorm())) * Eigen::Matrix<ctype, correctionSize, correctionSize>::Identity()
            - 4 * X * X.transpose();
      result.template bottomLeftCorner<1, correctionSize>() = 4 * X.transpose();

      // Upper hemisphere: adapt formulas so it is the stereographic projection from the south pole
      if (var[correctionSize] > 0) result.template bottomLeftCorner<1, correctionSize>() *= -1;

      // normalize the cols to make the orthogonal basis orthonormal
      result.colwise().normalize();

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
  bool operator==(const UnitVector<ctype2, d2> &v1, const UnitVector<ctype2, d2> &v2) {
    return v1.getValue() == v2.getValue();
  }

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

  template <typename ctype2, int d2>
  class RealTuple;

  template <typename ctype2, int d2, typename Scalar>
  requires std::is_arithmetic_v<Scalar>
  [[nodiscard]] RealTuple<ctype2, d2> operator*(const UnitVector<ctype2, d2> &rt, const Scalar &factor) {
    return RealTuple<ctype2, d2>(rt.getValue() * factor);
  }

  template <typename ctype2, int d2, typename Scalar>
  requires std::is_arithmetic_v<Scalar>
  [[nodiscard]] RealTuple<ctype2, d2> operator*(const Scalar &factor, const UnitVector<ctype2, d2> &rt) {
    return rt * factor;
  }

}  // namespace Ikarus
