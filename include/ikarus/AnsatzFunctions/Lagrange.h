/**
 * This implementation is taking from
 * and modified
 */

#ifndef IKARUS_LAGRANGE_H
#define IKARUS_LAGRANGE_H

#include <ikarus/utils/LinearAlgebraTypedefs.h>

#include <dune/common/power.hh>
/** \brief Lagrange shape functions of arbitrary order on the reference cube \f$[-1,1]^\text{dim}\f$
 *
 * This implementation is taken from
 https://gitlab.dune-project.org/core/dune-localfunctions/-/blob/master/dune/localfunctions/lagrange/lagrangecube.hh
 * and modified using Eigen Types

    Lagrange shape functions of arbitrary order have the property that
    \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

    \tparam ScalarType Type to represent the type of the coordinates
    \tparam dim Dimension of the domain cube
    \tparam k Polynomial order

    The ordering of the nodal shape functions are coordinate-wise and not counterclockwise or
 something similar!

  */
namespace Ikarus {

  template <class ScalarType, unsigned int dim, unsigned int k> class LagrangeCube {
  public:
    /** \brief Number of shape functions */
    static constexpr size_t sizeOfShapeFunctions = Dune::StaticPower<k + 1, dim>::power;

    /** \brief Type of the point where the ansatzfunction should be evaluated */
    using ParaMeterPointType = Ikarus::FixedVector<ScalarType, dim>;

    /** \brief VectorType of evaluated Ansatzfunctions */
    using VectorType = Ikarus::FixedVector<ScalarType, sizeOfShapeFunctions>;

    /** \brief Type of the Jacbian of the Ansatz Functions */
    using JacobianType = Ikarus::FixedMatrix<ScalarType, sizeOfShapeFunctions, dim>;

    static constexpr size_t size() { return sizeOfShapeFunctions; }

  private:
    // i-th Lagrange polynomial of degree k in one dimension
    static constexpr ScalarType getAnsatzFunctionImpl(unsigned int i, ScalarType x) {
      ScalarType result(1.0);
      for (unsigned int j = 0; j <= k; j++)
        if (j != i) result *= (k * x - j) / ((int)i - (int)j);
      return result;
    }

    // derivative of ith Lagrange polynomial of degree k in one dimension
    static ScalarType getAnsatzFunctionDerivativeImpl(unsigned int i, ScalarType x) {
      ScalarType result(0.0);

      for (unsigned int j = 0; j <= k; j++) {
        if (j != i) {
          ScalarType prod((k * 1.0) / ((int)i - (int)j));
          for (unsigned int l = 0; l <= k; l++)
            if (l != i && l != j) prod *= (k * x - l) / ((int)i - (int)l);
          result += prod;
        }
      }
      return result;
    }

    //! \brief transform x from -1..1 to 0..1
    static ParaMeterPointType transformPoint(const ParaMeterPointType& x) {
      return 0.5 * x.array() + 0.5;
    }

  public:
    //! \brief Evaluate all shape functions
    constexpr static VectorType evaluateFunction(const ParaMeterPointType& x01) {
      const ParaMeterPointType x(transformPoint(x01));
      VectorType N = VectorType::Zero();

      // Specialization for zero-order case
      if constexpr (k == 0)
        N[0] = 1;
      else if constexpr (k == 1)
        for (size_t i = 0; i < sizeOfShapeFunctions; i++) {
          N[i] = 1.0;
          for (unsigned int j = 0; j < dim; j++)
            // if j-th bit of i is set multiply with x[j], else with 1-x[j]
            N[i] *= (i & (1 << j)) ? x[j] : 1 - x[j];
        }
      else  // General case
        for (size_t i = 0; i < sizeOfShapeFunctions; i++) {
          // convert index i to multiindex
          std::array<unsigned int, dim> alpha(multiindex(i));

          // initialize product
          N[i] = 1.0;

          // dimension by dimension
          for (unsigned int j = 0; j < dim; j++) N[i] *= getAnsatzFunctionImpl(alpha[j], x[j]);
        }
      return N;  // transform back to -1..1
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference cube where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    static JacobianType evaluateJacobian(const ParaMeterPointType& x01) {
      const ParaMeterPointType x(transformPoint(x01));
      // Specialization for k==0
      if constexpr (k == 0) return JacobianType::Zero();
      JacobianType dN{};
      // Specialization for k==1
      if constexpr (k == 1) {
        // Loop over all shape functions
        for (size_t i = 0; i < sizeOfShapeFunctions; i++) {
          // Loop over all coordinate directions
          for (unsigned int j = 0; j < dim; j++) {
            // Initialize: the overall expression is a product
            // if j-th bit of i is set to 1, else -1
            dN(i, j) = (i & (1 << j)) ? 1 : -1;

            for (unsigned int l = 0; l < dim; l++) {
              if (j != l)
                // if l-th bit of i is set multiply with x[l], else with 1-x[l]
                dN(i, j) *= (i & (1 << l)) ? x[l] : 1 - x[l];
            }
          }
        }
      } else  // The general case
      {
        // Loop over all shape functions
        for (size_t i = 0; i < sizeOfShapeFunctions; i++) {
          // convert index i to multiindex
          std::array<unsigned int, dim> alpha(multiindex(i));

          // Loop over all coordinate directions
          for (unsigned int j = 0; j < dim; j++) {
            // Initialize: the overall expression is a product
            // if j-th bit of i is set to -1, else 1
            dN(i, j) = getAnsatzFunctionDerivativeImpl(alpha[j], x[j]);

            // rest of the product
            for (unsigned int l = 0; l < dim; l++)
              if (l != j) dN(i, j) *= getAnsatzFunctionImpl(alpha[l], x[l]);
          }
        }
      }
      return dN / 2;
    }

    static constexpr std::array<unsigned int, dim> multiindex(unsigned int i) {
      std::array<unsigned int, dim> alpha;
      for (unsigned int j = 0; j < dim; j++) {
        alpha[j] = i % (k + 1);
        i = i / (k + 1);
      }
      return alpha;
    }

    //! \brief Polynomial order of the shape functions
    static constexpr unsigned int order() { return k; }
  };

  //
  ///** \brief Lagrange shape functions of arbitrary order on the reference simplex
  // *
  // * This implementation is taken from
  // https://gitlab.dune-project.org/core/dune-localfunctions/-/blob/master/dune/localfunctions/lagrange/lagrangesimplex.hh
  // * and modified using eigen types
  //
  // Lagrange shape functions of arbitrary order have the property that
  //\f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.
  //
  //\tparam D Type to represent the field in the domain
  //\tparam R Type to represent the field in the range
  //\tparam dim Dimension of the domain simplex
  //\tparam k Polynomial order
  //*/
  // template<class ScalarType, unsigned int dim, unsigned int k, std::enable_if<!(dim!=3 &&
  // k>=2),bool> = true> class LagrangeSimplex
  //{
  // public:
  // public:
  //  /** \brief Number of shape functions */
  //  static constexpr size_t sizeOfShapeFunctions = Dune::binomial(k+dim,dim);
  //
  //  /** \brief Type of the point where the ansatzfunction should be evaluated */
  //  using ParaMeterPointType = Ikarus::FixedVector<ScalarType,dim>;
  //
  //  /** \brief VectorType of evaluated Ansatzfunctions */
  //  using VectorType = Ikarus::FixedVector<ScalarType,sizeOfShapeFunctions>;
  //
  //  /** \brief Type of the Jacbian of the Ansatz Functions */
  //  using JacobianType = Ikarus::FixedMatrix<ScalarType,sizeOfShapeFunctions,dim>;
  //
  //
  //  /** \brief Number of shape functions
  //   *
  //   * See https://en.wikipedia.org/wiki/Figurate_number for an explanation of the formula
  //   */
  //  static constexpr size_t size(){return sizeOfShapeFunctions;}
  //
  //  //! \brief Evaluate all shape functions
  //  constexpr static VectorType getAnsatzFunctions(const ParaMeterPointType& x)
  //  {
  //    VectorType N= VectorType::Zero();
  //    // Specialization for zero-order case
  //    if (k==0)
  //    {
  //      N[0] = 1;
  //      return N;
  //    }
  //
  //    // Specialization for first-order case
  //    if (k==1)
  //    {
  //      N[0] = 1.0;
  //      for (size_t i=0; i<dim; i++)
  //      {
  //        N[0]  -= x[i];
  //        N[i+1] = x[i];
  //      }
  //      return N;
  //    }
  //
  //    assert(k>=2);
  //
  //    auto lagrangeNode = [](unsigned int i) { return ((ScalarType)i)/k; };
  //
  //    if (dim==1)
  //    {
  //      for (unsigned int i=0; i<size(); i++)
  //      {
  //        N[i] = 1.0;
  //        for (unsigned int alpha=0; alpha<i; alpha++)
  //          N[i] *= (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
  //        for (unsigned int gamma=i+1; gamma<=k; gamma++)
  //          N[i] *= (x[0]-lagrangeNode(gamma))/(lagrangeNode(i)-lagrangeNode(gamma));
  //      }
  //      return;
  //    }
  //
  //    if (dim==2)
  //    {
  //      int n=0;
  //      for (unsigned int j=0; j<=k; j++)
  //        for (unsigned int i=0; i<=k-j; i++)
  //        {
  //          N[n] = 1.0;
  //          for (unsigned int alpha=0; alpha<i; alpha++)
  //            N[n] *= (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
  //          for (unsigned int beta=0; beta<j; beta++)
  //            N[n] *= (x[1]-lagrangeNode(beta))/(lagrangeNode(j)-lagrangeNode(beta));
  //          for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
  //            N[n] *=
  //            (lagrangeNode(gamma)-x[0]-x[1])/(lagrangeNode(gamma)-lagrangeNode(i)-lagrangeNode(j));
  //          n++;
  //        }
  //
  //      return;
  //    }
  //
  //    auto kx = x;
  //    kx *= k;
  //    unsigned int n = 0;
  //    unsigned int i[4];
  //    ScalarType factor[4];
  //    for (i[2] = 0; i[2] <= k; ++i[2])
  //    {
  //      factor[2] = 1.0;
  //      for (unsigned int j = 0; j < i[2]; ++j)
  //        factor[2] *= (kx[2]-j) / (i[2]-j);
  //      for (i[1] = 0; i[1] <= k - i[2]; ++i[1])
  //      {
  //        factor[1] = 1.0;
  //        for (unsigned int j = 0; j < i[1]; ++j)
  //          factor[1] *= (kx[1]-j) / (i[1]-j);
  //        for (i[0] = 0; i[0] <= k - i[1] - i[2]; ++i[0])
  //        {
  //          factor[0] = 1.0;
  //          for (unsigned int j = 0; j < i[0]; ++j)
  //            factor[0] *= (kx[0]-j) / (i[0]-j);
  //          i[3] = k - i[0] - i[1] - i[2];
  //          ScalarType kx3 = k - kx[0] - kx[1] - kx[2];
  //          factor[3] = 1.0;
  //          for (unsigned int j = 0; j < i[3]; ++j)
  //            factor[3] *= (kx3-j) / (i[3]-j);
  //          N[n++] = factor[0] * factor[1] * factor[2] * factor[3];
  //        }
  //      }
  //    }
  //    return N;
  //  }
  //
  //  /** \brief Evaluate Jacobian of all shape functions
  //   *
  //   * \param x Point in the reference simplex where to evaluation the Jacobians
  //   * \param[out] out The Jacobians of all shape functions at the point x
  //   */
  //  static JacobianType evaluateJacobian(const ParaMeterPointType& x)
  //  {
  //    JacobianType dN{};
  //
  //    // Specialization for k==0
  //    if (k==0)
  //    {
  //      return JacobianType::Zero();
  //    }
  //
  //    // Specialization for k==1
  //    if (k==1)
  //    {
  //      std::fill(out[0][0].begin(), out[0][0].end(), -1);
  //
  //      for (unsigned int i=0; i<dim; i++)
  //        for (unsigned int j=0; j<dim; j++)
  //          out[i+1][0][j] = (i==j);
  //
  //      return;
  //    }
  //
  //    auto lagrangeNode = [](unsigned int i) { return ((D)i)/k; };
  //
  //    // Specialization for dim==1
  //    if (dim==1)
  //    {
  //      for (unsigned int i=0; i<=k; i++)
  //      {
  //        // x_0 derivative
  //        out[i][0][0] = 0.0;
  //        R factor=1.0;
  //        for (unsigned int a=0; a<i; a++)
  //        {
  //          R product=factor;
  //          for (unsigned int alpha=0; alpha<i; alpha++)
  //            product *= (alpha==a) ? 1.0/(lagrangeNode(i)-lagrangeNode(alpha))
  //                                  :
  //                                  (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
  //          for (unsigned int gamma=i+1; gamma<=k; gamma++)
  //            product *= (lagrangeNode(gamma)-x[0])/(lagrangeNode(gamma)-lagrangeNode(i));
  //          out[i][0][0] += product;
  //        }
  //        for (unsigned int c=i+1; c<=k; c++)
  //        {
  //          R product=factor;
  //          for (unsigned int alpha=0; alpha<i; alpha++)
  //            product *= (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
  //          for (unsigned int gamma=i+1; gamma<=k; gamma++)
  //            product *= (gamma==c) ? -1.0/(lagrangeNode(gamma)-lagrangeNode(i))
  //                                  :
  //                                  (lagrangeNode(gamma)-x[0])/(lagrangeNode(gamma)-lagrangeNode(i));
  //          out[i][0][0] += product;
  //        }
  //      }
  //      return;
  //    }
  //
  //    if (dim==2)
  //    {
  //      int n=0;
  //      for (unsigned int j=0; j<=k; j++)
  //        for (unsigned int i=0; i<=k-j; i++)
  //        {
  //          // x_0 derivative
  //          out[n][0][0] = 0.0;
  //          R factor=1.0;
  //          for (unsigned int beta=0; beta<j; beta++)
  //            factor *= (x[1]-lagrangeNode(beta))/(lagrangeNode(j)-lagrangeNode(beta));
  //          for (unsigned int a=0; a<i; a++)
  //          {
  //            R product=factor;
  //            for (unsigned int alpha=0; alpha<i; alpha++)
  //              if (alpha==a)
  //                product *= D(1)/(lagrangeNode(i)-lagrangeNode(alpha));
  //              else
  //                product *= (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
  //            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
  //              product *=
  //              (lagrangeNode(gamma)-x[0]-x[1])/(lagrangeNode(gamma)-lagrangeNode(i)-lagrangeNode(j));
  //            out[n][0][0] += product;
  //          }
  //          for (unsigned int c=i+j+1; c<=k; c++)
  //          {
  //            R product=factor;
  //            for (unsigned int alpha=0; alpha<i; alpha++)
  //              product *= (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
  //            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
  //              if (gamma==c)
  //                product *= -D(1)/(lagrangeNode(gamma)-lagrangeNode(i)-lagrangeNode(j));
  //              else
  //                product *=
  //                (lagrangeNode(gamma)-x[0]-x[1])/(lagrangeNode(gamma)-lagrangeNode(i)-lagrangeNode(j));
  //            out[n][0][0] += product;
  //          }
  //
  //          // x_1 derivative
  //          out[n][0][1] = 0.0;
  //          factor = 1.0;
  //          for (unsigned int alpha=0; alpha<i; alpha++)
  //            factor *= (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
  //          for (unsigned int b=0; b<j; b++)
  //          {
  //            R product=factor;
  //            for (unsigned int beta=0; beta<j; beta++)
  //              if (beta==b)
  //                product *= D(1)/(lagrangeNode(j)-lagrangeNode(beta));
  //              else
  //                product *= (x[1]-lagrangeNode(beta))/(lagrangeNode(j)-lagrangeNode(beta));
  //            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
  //              product *=
  //              (lagrangeNode(gamma)-x[0]-x[1])/(lagrangeNode(gamma)-lagrangeNode(i)-lagrangeNode(j));
  //            out[n][0][1] += product;
  //          }
  //          for (unsigned int c=i+j+1; c<=k; c++)
  //          {
  //            R product=factor;
  //            for (unsigned int beta=0; beta<j; beta++)
  //              product *= (x[1]-lagrangeNode(beta))/(lagrangeNode(j)-lagrangeNode(beta));
  //            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
  //              if (gamma==c)
  //                product *= -D(1)/(lagrangeNode(gamma)-lagrangeNode(i)-lagrangeNode(j));
  //              else
  //                product *=
  //                (lagrangeNode(gamma)-x[0]-x[1])/(lagrangeNode(gamma)-lagrangeNode(i)-lagrangeNode(j));
  //            out[n][0][1] += product;
  //          }
  //
  //          n++;
  //        }
  //
  //      return;
  //    }
  //
  //    if (dim!=3)
  //      DUNE_THROW(NotImplemented, "LagrangeSimplexLocalBasis only implemented for dim==3!");
  //
  //    // Specialization for arbitrary order and dim==3
  //    typename Traits::DomainType kx = x;
  //    kx *= k;
  //    unsigned int n = 0;
  //    unsigned int i[4];
  //    R factor[4];
  //    for (i[2] = 0; i[2] <= k; ++i[2])
  //    {
  //      factor[2] = 1.0;
  //      for (unsigned int j = 0; j < i[2]; ++j)
  //        factor[2] *= (kx[2]-j) / (i[2]-j);
  //      for (i[1] = 0; i[1] <= k - i[2]; ++i[1])
  //      {
  //        factor[1] = 1.0;
  //        for (unsigned int j = 0; j < i[1]; ++j)
  //          factor[1] *= (kx[1]-j) / (i[1]-j);
  //        for (i[0] = 0; i[0] <= k - i[1] - i[2]; ++i[0])
  //        {
  //          factor[0] = 1.0;
  //          for (unsigned int j = 0; j < i[0]; ++j)
  //            factor[0] *= (kx[0]-j) / (i[0]-j);
  //          i[3] = k - i[0] - i[1] - i[2];
  //          D kx3 = k - kx[0] - kx[1] - kx[2];
  //          R sum3 = 0.0;
  //          factor[3] = 1.0;
  //          for (unsigned int j = 0; j < i[3]; ++j)
  //            factor[3] /= i[3] - j;
  //          R prod_all = factor[0] * factor[1] * factor[2] * factor[3];
  //          for (unsigned int j = 0; j < i[3]; ++j)
  //          {
  //            R prod = prod_all;
  //            for (unsigned int l = 0; l < i[3]; ++l)
  //              if (j == l)
  //                prod *= -R(k);
  //              else
  //                prod *= kx3 - l;
  //            sum3 += prod;
  //          }
  //          for (unsigned int j = 0; j < i[3]; ++j)
  //            factor[3] *= kx3 - j;
  //          for (unsigned int m = 0; m < 3; ++m)
  //          {
  //            out[n][0][m] = sum3;
  //            for (unsigned int j = 0; j < i[m]; ++j)
  //            {
  //              R prod = factor[3];
  //              for (unsigned int p = 0; p < 3; ++p)
  //              {
  //                if (m == p)
  //                  for (unsigned int l = 0; l < i[p]; ++l)
  //                    prod *= (j==l) ? R(k) / (i[p]-l) : (kx[p]-l) / (i[p]-l);
  //                else
  //                  prod *= factor[p];
  //              }
  //              out[n][0][m] += prod;
  //            }
  //          }
  //          n++;
  //        }
  //      }
  //    }
  //  }
  //
  //  /** \brief Evaluate partial derivatives of any order of all shape functions
  //   *
  //   * \param order Order of the partial derivatives, in the classic multi-index notation
  //   * \param in Position where to evaluate the derivatives
  //   * \param[out] out The desired partial derivatives
  //   */
  //  void partial(const std::array<unsigned int,dim>& order,
  //               const typename Traits::DomainType& in,
  //               std::vector<typename Traits::RangeType>& out) const
  //  {
  //    auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
  //
  //    out.resize(size());
  //
  //    if (totalOrder == 0) {
  //      evaluateFunction(in, out);
  //      return;
  //    }
  //
  //    if (k==0)
  //    {
  //      out[0] = 0;
  //      return;
  //    }
  //
  //    if (k==1)
  //    {
  //      if (totalOrder==1)
  //      {
  //        auto direction = std::find(order.begin(), order.end(), 1);
  //
  //        out[0] = -1;
  //        for (unsigned int i=0; i<dim; i++)
  //          out[i+1] = (i==(direction-order.begin()));
  //      }
  //      else  // all higher order derivatives are zero
  //        std::fill(out.begin(), out.end(), 0);
  //      return;
  //    }
  //
  //    if (dim==2)
  //    {
  //      auto lagrangeNode = [](unsigned int i) { return ((D)i)/k; };
  //
  //      // Helper method: Return a single Lagrangian factor of l_ij evaluated at x
  //      auto lagrangianFactor = [&lagrangeNode]
  //          (const int no, const int i, const int j, const typename Traits::DomainType& x)
  //          -> typename Traits::RangeType
  //      {
  //        if ( no < i)
  //          return (x[0]-lagrangeNode(no))/(lagrangeNode(i)-lagrangeNode(no));
  //        if (no < i+j)
  //          return (x[1]-lagrangeNode(no-i))/(lagrangeNode(j)-lagrangeNode(no-i));
  //        return
  //        (lagrangeNode(no+1)-x[0]-x[1])/(lagrangeNode(no+1)-lagrangeNode(i)-lagrangeNode(j));
  //      };
  //
  //      // Helper method: Return the derivative of a single Lagrangian factor of l_ij evaluated at
  //      x
  //      // direction: Derive in x-direction if this is 0, otherwise derive in y direction
  //      auto lagrangianFactorDerivative = [&lagrangeNode]
  //          (const int direction, const int no, const int i, const int j, const typename
  //          Traits::DomainType& x)
  //          -> typename Traits::RangeType
  //      {
  //        using T = typename Traits::RangeType;
  //        if ( no < i)
  //          return (direction == 0) ? T(1.0/(lagrangeNode(i)-lagrangeNode(no))) : T(0);
  //
  //        if (no < i+j)
  //          return (direction == 0) ? T(0) : T(1.0/(lagrangeNode(j)-lagrangeNode(no-i)));
  //
  //        return -1.0/(lagrangeNode(no+1)-lagrangeNode(i)-lagrangeNode(j));
  //      };
  //
  //      if (totalOrder==1)
  //      {
  //        int direction = std::find(order.begin(), order.end(), 1)-order.begin();
  //
  //        int n=0;
  //        for (unsigned int j=0; j<=k; j++)
  //        {
  //          for (unsigned int i=0; i<=k-j; i++, n++)
  //          {
  //            out[n] = 0.0;
  //            for (unsigned int no1=0; no1 < k; no1++)
  //            {
  //              R factor = lagrangianFactorDerivative(direction, no1, i, j, in);
  //              for (unsigned int no2=0; no2 < k; no2++)
  //                if (no1 != no2)
  //                  factor *= lagrangianFactor(no2, i, j, in);
  //
  //              out[n] += factor;
  //            }
  //          }
  //        }
  //        return;
  //      }
  //
  //      if (totalOrder==2)
  //      {
  //        std::array<int,2> directions;
  //        unsigned int counter = 0;
  //        auto nonconstOrder = order;  // need a copy that I can modify
  //        for (int i=0; i<2; i++)
  //        {
  //          while (nonconstOrder[i])
  //          {
  //            directions[counter++] = i;
  //            nonconstOrder[i]--;
  //          }
  //        }
  //
  //        //f = prod_{i} f_i -> dxa dxb f = sum_{i} {dxa f_i sum_{k \neq i} dxb f_k prod_{l \neq
  //        k,i} f_l int n=0; for (unsigned int j=0; j<=k; j++)
  //        {
  //          for (unsigned int i=0; i<=k-j; i++, n++)
  //          {
  //            R res = 0.0;
  //
  //            for (unsigned int no1=0; no1 < k; no1++)
  //            {
  //              R factor1 = lagrangianFactorDerivative(directions[0], no1, i, j, in);
  //              for (unsigned int no2=0; no2 < k; no2++)
  //              {
  //                if (no1 == no2)
  //                  continue;
  //                R factor2 = factor1*lagrangianFactorDerivative(directions[1], no2, i, j, in);
  //                for (unsigned int no3=0; no3 < k; no3++)
  //                {
  //                  if (no3 == no1 || no3 == no2)
  //                    continue;
  //                  factor2 *= lagrangianFactor(no3, i, j, in);
  //                }
  //                res += factor2;
  //              }
  //            }
  //            out[n] = res;
  //          }
  //        }
  //
  //        return;
  //      }  // totalOrder==2
  //
  //    }   // dim==2
  //
  //    DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
  //  }
  //
  //  //! \brief Polynomial order of the shape functions
  //  static constexpr unsigned int order ()
  //  {
  //    return k;
  //  }
  //};

}  // namespace Ikarus

#endif  // IKARUS_LAGRANGE_H
