//
// Created by lex on 04/02/2022.
//

#pragma once
#include <Eigen/Core>
namespace Ikarus
{
  template<typename ST, int size> requires (size>0 and size <= 3)
  auto toVoigt(const Eigen::Matrix<ST,size,size>& E)
  {
    Eigen::Vector<ST, (size * (size + 1)) / 2> EVoigt;
    EVoigt.setZero();
    for (int i = 0; i < size; ++i)
      EVoigt(i) = E(i, i);

    if constexpr (size > 1) EVoigt(size) = E(0, 1) * 2;
    if constexpr (size > 2) {
      EVoigt(size + 1) = E(0, 2) * 2;
      EVoigt(size + 2) = E(1, 2) * 2;
    }
   return EVoigt;
  }

  template <typename LocalView>
  struct TraitsFromLocalView {
    using GridEntity = typename LocalView::Element;
    /** \brief Dimension of the world space */
    static constexpr int worlddim = GridEntity::Geometry::coorddimension;

    /** \brief Dimension of the geometry */
    static constexpr int mydim = GridEntity::mydimension;

    /** \brief Dimension of the grid */
    static constexpr int dimension = GridEntity::dimension;

    /** \brief Type of the internal forces */
    using VectorType = Eigen::VectorXd;

    /** \brief Type of the stiffness matrix */
    using MatrixType = Eigen::MatrixXd;

    /** \brief Type of the stiffness matrix */
    using ScalarType = typename GridEntity::Geometry::ctype;
  };

  struct EmodulAndNu {
    double emodul;
    double nu;
  };

  template <typename ValuePair>
  struct ConvertLameConstants {
    ConvertLameConstants(ValuePair&& p_vp) : vp(p_vp) {}


    double toLambda() {
      if constexpr (std::is_same_v<ValuePair, EmodulAndNu>) {
        auto& E  = vp.emodul;
        auto& nu = vp.nu;
        return E * nu / ((1 + nu) * (1 - 2 * nu));
      }
      else assert(false && "Your LameParameter request is not implemented");
    }
    ValuePair vp;
  };

  //  template <typename ValuePair, typename ConvertTo>
  //  double convertLameConstants(ValuePair&& p_vp, ConvertTo&& convertTo )
  //  {
  //    return ConvertLameConstants(p_vp).t
  //  }


}


