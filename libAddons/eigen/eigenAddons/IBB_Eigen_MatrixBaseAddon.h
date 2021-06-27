//
// Created by alex on 05.05.20.
//

#ifndef NUMPRO_IBB_EIGEN_MATRIXBASEADDON_H
#define NUMPRO_IBB_EIGEN_MATRIXBASEADDON_H

friend void PrintTo(const Derived& m, ::std::ostream* o) { *o << "\n" << m; }

/** Berechnet die Inverse eines 3x3-Array2D
 *
 * \author Alex Müller
 * \version 22.01.2020
 */
inline Derived inverse_33() { return this->inverse(); }

/** partially copied from\src\ThirdParty\eigen3\Eigen\src\Geometry\OrthoMethods.h but removed assert
 * (ONLY FOR LEGACY DO NOT USE THIS)
 *
 * \author Alex Müller
 * \version 09.05.2020
 */
template <typename OtherDerived>
Derived cross_product(const MatrixBase<OtherDerived>& other) const {
  //  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived,3)
  //  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(OtherDerived,3)

  typedef typename Derived::Scalar Scalar;

  typename internal::nested_eval<Derived, 2>::type lhs(derived());
  typename internal::nested_eval<OtherDerived, 2>::type rhs(other.derived());

  Matrix<Scalar, 3, 1> thisSurrogate = (Matrix<Scalar, 3, 1>() << lhs.coeff(0), lhs.coeff(1), lhs.coeff(2)).finished();
  Matrix<Scalar, 3, 1> otherSurrogate = (Matrix<Scalar, 3, 1>() << rhs.coeff(0), rhs.coeff(1), rhs.coeff(2)).finished();
  return thisSurrogate.cross(otherSurrogate);
}

/** Berechnet die Loesung eines LGS mit einem Array2D
 *
 * \author Alex Müller
 * \version 22.01.2020
 */
template <typename OtherDerived>
inline OtherDerived solve(const MatrixBase<OtherDerived>& rhs) const {
  return this->partialPivLu().solve(rhs);
}

/** nicht implementiert!!
 *
 * \author Malte von Scheven
 * \version 25.08.2011
 */
void print_mask() {
  //  NP_log("\n");  //TODO Alex/Tobi: Funktion implementieren oder rauswerfen(leer lassen)?
  //  for (int z=0; z<this->cols(); z++)
  //  {
  //    for (int s=0; s<this->rows(); s++)
  //      NP_logf("%1s", fabs(this(z,s) ) > 1e-8? "X":" " );
  //    NP_log("\n");
  //  }
  //  NP_log("\n");
}

// Derived orthonormalizeMatrixColumns() {
////Gram Schmidt Ortho
//  Derived Q = *this;
//
//  Q.col(0).normalize();
//
//  for (int colIndex = 1; colIndex < Q.cols(); colIndex++) {
//    Q.col(colIndex) -= Q.leftCols(colIndex) * (Q.leftCols(colIndex).transpose() *
//    (*this).col(colIndex)); Q.col(colIndex).normalize();
//  }
//
//  return Q;
//}

#endif  // NUMPRO_IBB_EIGEN_MATRIXBASEADDON_H
