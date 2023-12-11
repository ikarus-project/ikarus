// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#if HAVE_DUNE_LOCALFEFUNCTIONS
#  include <dune/localfefunctions/derivativetransformators.hh>
#  include <dune/localfefunctions/meta.hh>

#  include <ikarus/finiteelements/ferequirements.hh>
#  include <ikarus/utils/eigendunetransformations.hh>

namespace Ikarus {

  template <typename Geometry>
  Eigen::Matrix3d calcTransformationMatrix2D(const Geometry& geometry) {
    const auto& referenceElement = Dune::ReferenceElements<double, 2>::general(geometry.type());
    const auto quadPos0          = referenceElement.position(0, 0);

    const auto jacobianinvT0 = toEigen(geometry.jacobianInverseTransposed(quadPos0));
    const auto detJ0         = geometry.integrationElement(quadPos0);

    auto jaco = (jacobianinvT0).inverse().eval();
    auto J11  = jaco(0, 0);
    auto J12  = jaco(0, 1);
    auto J21  = jaco(1, 0);
    auto J22  = jaco(1, 1);

    Eigen::Matrix3d T0;
    // clang-format off
    T0 <<      J11 * J11, J12 * J12,                   J11 * J12,
               J21 * J21, J22 * J22,                   J21 * J22,
         2.0 * J11 * J21, 2.0 * J12 * J22, J21 * J12 + J11 * J22;
    // clang-format on
    return T0.inverse() * detJ0;
  }

  template <typename Geometry>
  Eigen::Matrix<double, 6, 6> calcTransformationMatrix3D(const Geometry& geometry) {
    const auto& referenceElement = Dune::ReferenceElements<double, 3>::general(geometry.type());
    const auto quadPos0          = referenceElement.position(0, 0);

    const auto jacobianinvT0 = toEigen(geometry.jacobianInverseTransposed(quadPos0));
    const auto detJ0         = geometry.integrationElement(quadPos0);

    auto jaco = (jacobianinvT0).inverse().eval();
    auto J11  = jaco(0, 0);
    auto J12  = jaco(0, 1);
    auto J13  = jaco(0, 2);
    auto J21  = jaco(1, 0);
    auto J22  = jaco(1, 1);
    auto J23  = jaco(1, 2);
    auto J31  = jaco(2, 0);
    auto J32  = jaco(2, 1);
    auto J33  = jaco(2, 2);

    Eigen::Matrix<double, 6, 6> T0;
    // clang-format off
    T0 <<      J11 * J11,       J12 * J12,       J13 * J13,             J11 * J12,             J11 * J13,             J12 * J13,
               J21 * J21,       J22 * J22,       J23 * J23,             J21 * J22,             J21 * J23,             J22 * J23,
               J31 * J31,       J32 * J32,       J33 * J33,             J31 * J32,             J31 * J33,             J32 * J33,
         2.0 * J11 * J21, 2.0 * J12 * J22, 2.0 * J13 * J23, J11 * J22 + J21 * J12, J11 * J23 + J21 * J13, J12 * J23 + J22 * J13,
         2.0 * J11 * J31, 2.0 * J12 * J32, 2.0 * J13 * J33, J11 * J32 + J31 * J12, J11 * J33 + J31 * J13, J12 * J33 + J32 * J13,
         2.0 * J31 * J21, 2.0 * J32 * J22, 2.0 * J33 * J23, J31 * J22 + J21 * J32, J31 * J23 + J21 * J33, J32 * J23 + J22 * J33;
    // clang-format on

    return T0.inverse() * detJ0;
  }

  template <typename Geometry>
  struct EASQ1E4 {
    static constexpr int strainSize         = 3;
    static constexpr int enhancedStrainSize = 4;
    using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;

    EASQ1E4() = default;
    explicit EASQ1E4(const Geometry& geometry_)
        : geometry{std::make_shared<Geometry>(geometry_)},
          T0InverseTransformed{calcTransformationMatrix2D(geometry_)} {}

    auto calcM(const Dune::FieldVector<double, 2>& quadPos) const {
      MType M;
      M.setZero(strainSize, enhancedStrainSize);
      const double xi   = quadPos[0];
      const double eta  = quadPos[1];
      M(0, 0)           = 2 * xi - 1.0;
      M(1, 1)           = 2 * eta - 1.0;
      M(2, 2)           = 2 * xi - 1.0;
      M(2, 3)           = 2 * eta - 1.0;
      const double detJ = geometry->integrationElement(quadPos);
      M                 = T0InverseTransformed / detJ * M;
      return M;
    }

    std::shared_ptr<Geometry> geometry;
    Eigen::Matrix3d T0InverseTransformed;
  };

  template <typename Geometry>
  struct EASQ1E5 {
    static constexpr int strainSize         = 3;
    static constexpr int enhancedStrainSize = 5;
    using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;

    EASQ1E5() = default;
    explicit EASQ1E5(const Geometry& geometry_)
        : geometry{std::make_shared<Geometry>(geometry_)},
          T0InverseTransformed{calcTransformationMatrix2D(geometry_)} {}

    auto calcM(const Dune::FieldVector<double, 2>& quadPos) const {
      MType M;
      M.setZero();
      const double xi   = quadPos[0];
      const double eta  = quadPos[1];
      M(0, 0)           = 2 * xi - 1.0;
      M(1, 1)           = 2 * eta - 1.0;
      M(2, 2)           = 2 * xi - 1.0;
      M(2, 3)           = 2 * eta - 1.0;
      M(2, 4)           = (2 * xi - 1.0) * (2 * eta - 1.0);
      const double detJ = geometry->integrationElement(quadPos);
      M                 = T0InverseTransformed / detJ * M;
      return M;
    }

    std::shared_ptr<Geometry> geometry;
    Eigen::Matrix3d T0InverseTransformed;
  };

  template <typename Geometry>
  struct EASQ1E7 {
    static constexpr int strainSize         = 3;
    static constexpr int enhancedStrainSize = 7;
    using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;

    EASQ1E7() = default;
    explicit EASQ1E7(const Geometry& geometry_)
        : geometry{std::make_shared<Geometry>(geometry_)},
          T0InverseTransformed{calcTransformationMatrix2D(geometry_)} {}

    auto calcM(const Dune::FieldVector<double, 2>& quadPos) const {
      MType M;
      M.setZero();
      const double xi   = quadPos[0];
      const double eta  = quadPos[1];
      M(0, 0)           = 2 * xi - 1.0;
      M(1, 1)           = 2 * eta - 1.0;
      M(2, 2)           = 2 * xi - 1.0;
      M(2, 3)           = 2 * eta - 1.0;
      M(0, 4)           = (2 * xi - 1.0) * (2 * eta - 1.0);
      M(1, 5)           = (2 * xi - 1.0) * (2 * eta - 1.0);
      M(2, 6)           = (2 * xi - 1.0) * (2 * eta - 1.0);
      const double detJ = geometry->integrationElement(quadPos);
      M                 = T0InverseTransformed / detJ * M;
      return M;
    }

    std::shared_ptr<Geometry> geometry;
    Eigen::Matrix3d T0InverseTransformed;
  };

  template <typename Geometry>
  struct EASH1E9 {
    static constexpr int strainSize         = 6;
    static constexpr int enhancedStrainSize = 9;
    using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;

    EASH1E9() = default;
    explicit EASH1E9(const Geometry& geometry_)
        : geometry{std::make_shared<Geometry>(geometry_)},
          T0InverseTransformed{calcTransformationMatrix3D(geometry_)} {}

    auto calcM(const Dune::FieldVector<double, 3>& quadPos) const {
      MType M;
      M.setZero();
      const double xi   = quadPos[0];
      const double eta  = quadPos[1];
      const double zeta = quadPos[2];
      M(0, 0)           = 2 * xi - 1.0;
      M(1, 1)           = 2 * eta - 1.0;
      M(2, 2)           = 2 * zeta - 1.0;
      M(3, 3)           = 2 * xi - 1.0;
      M(3, 4)           = 2 * eta - 1.0;
      M(4, 5)           = 2 * xi - 1.0;
      M(4, 6)           = 2 * zeta - 1.0;
      M(5, 7)           = 2 * eta - 1.0;
      M(5, 8)           = 2 * zeta - 1.0;
      const double detJ = geometry->integrationElement(quadPos);
      M                 = T0InverseTransformed / detJ * M;
      return M;
    }
    std::shared_ptr<Geometry> geometry;
    Eigen::Matrix<double, 6, 6> T0InverseTransformed;
  };

  template <typename Geometry>
  struct EASH1E21 {
    static constexpr int strainSize         = 6;
    static constexpr int enhancedStrainSize = 21;
    using MType                             = Eigen::Matrix<double, strainSize, enhancedStrainSize>;

    EASH1E21() = default;
    explicit EASH1E21(const Geometry& geometry_)
        : geometry{std::make_shared<Geometry>(geometry_)},
          T0InverseTransformed{calcTransformationMatrix3D(geometry_)} {}

    auto calcM(const Dune::FieldVector<double, 3>& quadPos) const {
      MType M;
      M.setZero();
      const double xi   = quadPos[0];
      const double eta  = quadPos[1];
      const double zeta = quadPos[2];
      M(0, 0)           = 2 * xi - 1.0;
      M(1, 1)           = 2 * eta - 1.0;
      M(2, 2)           = 2 * zeta - 1.0;
      M(3, 3)           = 2 * xi - 1.0;
      M(3, 4)           = 2 * eta - 1.0;
      M(4, 5)           = 2 * xi - 1.0;
      M(4, 6)           = 2 * zeta - 1.0;
      M(5, 7)           = 2 * eta - 1.0;
      M(5, 8)           = 2 * zeta - 1.0;

      M(3, 9)  = (2 * xi - 1.0) * (2 * zeta - 1.0);
      M(3, 10) = (2 * eta - 1.0) * (2 * zeta - 1.0);
      M(4, 11) = (2 * xi - 1.0) * (2 * eta - 1.0);
      M(4, 12) = (2 * eta - 1.0) * (2 * zeta - 1.0);
      M(5, 13) = (2 * xi - 1.0) * (2 * eta - 1.0);
      M(5, 14) = (2 * xi - 1.0) * (2 * zeta - 1.0);

      M(0, 15) = (2 * xi - 1.0) * (2 * eta - 1.0);
      M(0, 16) = (2 * xi - 1.0) * (2 * zeta - 1.0);
      M(1, 17) = (2 * xi - 1.0) * (2 * eta - 1.0);
      M(1, 18) = (2 * eta - 1.0) * (2 * zeta - 1.0);
      M(2, 19) = (2 * xi - 1.0) * (2 * zeta - 1.0);
      M(2, 20) = (2 * eta - 1.0) * (2 * zeta - 1.0);

      const double detJ = geometry->integrationElement(quadPos);
      M                 = T0InverseTransformed / detJ * M;
      return M;
    }

    std::shared_ptr<Geometry> geometry;
    Eigen::Matrix<double, 6, 6> T0InverseTransformed;
  };

  template <typename Geometry>
  using EAS2DVariant = std::variant<std::monostate, EASQ1E4<Geometry>, EASQ1E5<Geometry>, EASQ1E7<Geometry>>;
  template <typename Geometry>
  using EAS3DVariant = std::variant<std::monostate, EASH1E9<Geometry>, EASH1E21<Geometry>>;

  template <typename DisplacementBasedElement>
  class EnhancedAssumedStrains : public DisplacementBasedElement {
  public:
    using FERequirementType      = typename DisplacementBasedElement::FERequirementType;
    using ResultRequirementsType = typename DisplacementBasedElement::ResultRequirementsType;
    using LocalView              = typename DisplacementBasedElement::LocalView;
    using GridView               = typename DisplacementBasedElement::GridView;
    using Traits                 = typename DisplacementBasedElement::Traits;
    using DisplacementBasedElement::localView;

    // Disabling this forwarding constructor if the argument provided is EnhancedAssumedStrains itself, to forward the
    // calls to the implicit copy constructor
    template <typename... Args>
    requires(not std::is_same_v<std::remove_cvref_t<std::tuple_element_t<0, std::tuple<Args...>>>,
                                EnhancedAssumedStrains>) explicit EnhancedAssumedStrains(Args&&... args)
        : DisplacementBasedElement(std::forward<Args>(args)...) {}

    inline double calculateScalar(const FERequirementType& par) const {
      if (isDisplacementBased()) return DisplacementBasedElement::calculateScalar(par);
      DUNE_THROW(Dune::NotImplemented,
                 "EAS element do not support any scalar calculations, i.e. they are not derivable from a potential");
    }

    bool isDisplacementBased() const { return std::holds_alternative<std::monostate>(easVariant_); }

    inline void calculateVector(const FERequirementType& par, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(par, force);
    }

    const auto& easVariant() const { return easVariant_; }
    auto getNumberOfEASParameters() const {
      return std::visit(
          [&]<typename EAST>(const EAST&) {
            if constexpr (std::is_same_v<std::monostate, EAST>)
              return 0;
            else
              return EAST::enhancedStrainSize;
          },
          easVariant_);
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::template MatrixType<> K) const {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      /// fill K with displacement-based stiffness.
      /// It is assumed to be assembled block-wise on element level.
      /// This means the displacements x,y,z of node I are grouped together.
      DisplacementBasedElement::calculateMatrix(par, K);

      if (isDisplacementBased()) return;

      std::visit(
          [&]<typename EAST>(const EAST& easFunction) {
            if constexpr (not std::is_same_v<std::monostate, EAST>) {
              constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
              Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize> D;
              calculateDAndLMatrix(easFunction, par, D, L);

              K.template triangularView<Eigen::Upper>() -= L.transpose() * D.inverse() * L;
              K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
            }
          },
          easVariant_);
    }

    void calculateAt(const ResultRequirementsType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::Indices;
      using namespace Dune::DerivativeDirections;
      using namespace Dune;

      DisplacementBasedElement::calculateAt(req, local, result);

      if (isDisplacementBased()) return;

      const auto& d        = req.getGlobalSolution(Ikarus::FESolutions::displacement);
      const auto C         = DisplacementBasedElement::getMaterialTangentFunction(req.getFERequirements());
      const auto& numNodes = DisplacementBasedElement::numberOfNodes;

      Eigen::VectorXd disp(localView().size());
      for (auto i = 0U; i < numNodes; ++i)
        for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
          disp[i * Traits::mydim + k2] = d[localView().index(localView().tree().child(k2).localIndex(i))[0]];

      std::visit(
          [&]<typename EAST>(const EAST& easFunction) {
            if constexpr (not std::is_same_v<std::monostate, EAST>) {
              constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
              Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize> D;
              calculateDAndLMatrix(easFunction, req.getFERequirements(), D, L);
              const auto alpha = (-D.inverse() * L * disp).eval();
              const auto M     = easFunction.calcM(local);
              const auto CEval = C(local);
              auto easStress   = (CEval * M * alpha).eval();
              typename ResultTypeMap<double>::ResultArray resultVector;
              if (req.isResultRequested(ResultType::linearStress)) {
                resultVector.resize(3, 1);
                resultVector = result.getResult(ResultType::linearStress) + easStress;
                result.insertOrAssignResult(ResultType::linearStress, resultVector);
              } else
                DUNE_THROW(Dune::NotImplemented, "The requested result type is NOT implemented.");
            }
          },
          easVariant_);
    }

    void setEASType(int numberOfEASParameters_) {
      const auto& numNodes = DisplacementBasedElement::numberOfNodes;
      if (not((numNodes == 4 and Traits::mydim == 2) or (numNodes == 8 and Traits::mydim == 3))
          and (not isDisplacementBased()))
        DUNE_THROW(Dune::NotImplemented, "EAS only supported for Q1 or H1 elements");

      if constexpr (Traits::mydim == 2) {
        switch (numberOfEASParameters_) {
          case 0:
            easVariant_ = std::monostate();
            break;
          case 4:
            easVariant_ = EASQ1E4(localView().element().geometry());
            break;
          case 5:
            easVariant_ = EASQ1E5(localView().element().geometry());
            break;
          case 7:
            easVariant_ = EASQ1E7(localView().element().geometry());
            break;
          default:
            DUNE_THROW(Dune::NotImplemented, "The given EAS parameters are not available for the 2D case.");
            break;
        }
      } else if constexpr (Traits::mydim == 3) {
        switch (numberOfEASParameters_) {
          case 0:
            easVariant_ = std::monostate();
            break;
          case 9:
            easVariant_ = EASH1E9(localView().element().geometry());
            break;
          case 21:
            easVariant_ = EASH1E21(localView().element().geometry());
            break;
          default:
            DUNE_THROW(Dune::NotImplemented, "The given EAS parameters are not available for the 3D case.");
            break;
        }
      }
    }

  protected:
    template <typename ScalarType>
    inline ScalarType calculateScalarImpl(const FERequirementType& par,
                                          const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                          = std::nullopt) const {
      if (isDisplacementBased()) return DisplacementBasedElement::template calculateScalarImpl<ScalarType>(par, dx);
      DUNE_THROW(Dune::NotImplemented,
                 "EAS element do not support any scalar calculations, i.e. they are not derivable from a potential");
    }

    template <typename ScalarType>
    void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      DisplacementBasedElement::calculateVectorImpl(par, force, dx);
      if (isDisplacementBased()) return;
      using namespace Dune;
      const auto& d       = par.getGlobalSolution(Ikarus::FESolutions::displacement);
      auto strainFunction = DisplacementBasedElement::getStrainFunction(par, dx);
      Eigen::VectorX<ScalarType> disp(localView().size());
      const auto& numNodes = DisplacementBasedElement::numberOfNodes;

      // FIXME this should not be needed in the future strainFunction should be able to hand out this vector
      if (dx)
        for (auto i = 0U; i < numNodes; ++i)
          for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
            disp[i * Traits::mydim + k2] = dx.value()[i * Traits::mydim + k2]
                                           + d[localView().index(localView().tree().child(k2).localIndex(i))[0]];
      else
        for (auto i = 0U; i < numNodes; ++i)
          for (auto k2 = 0U; k2 < Traits::mydim; ++k2)
            disp[i * Traits::mydim + k2] = d[localView().index(localView().tree().child(k2).localIndex(i))[0]];

      using namespace Dune::DerivativeDirections;

      auto C         = DisplacementBasedElement::getMaterialTangentFunction(par);
      const auto geo = localView().element().geometry();

      // Internal forces from enhanced strains
      std::visit(
          [&]<typename EAST>(const EAST& easFunction) {
            if constexpr (not std::is_same_v<std::monostate, EAST>) {
              constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
              Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize> D;
              calculateDAndLMatrix(easFunction, par, D, L);

              const auto alpha = (-D.inverse() * L * disp).eval();

              for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
                const auto M            = easFunction.calcM(gp.position());
                const double intElement = geo.integrationElement(gp.position()) * gp.weight();
                const auto CEval        = C(gpIndex);
                auto stresses           = (CEval * M * alpha).eval();
                for (size_t i = 0; i < numNodes; ++i) {
                  const auto bopI = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
                  force.template segment<Traits::mydim>(Traits::mydim * i) += bopI.transpose() * stresses * intElement;
                }
              }
            }
          },
          easVariant_);
    }

  private:
    using EAS2DVariantImpl = EAS2DVariant<typename LocalView::Element::Geometry>;
    using EAS3DVariantImpl = EAS3DVariant<typename LocalView::Element::Geometry>;
    std::conditional_t<Traits::mydim == 2, EAS2DVariantImpl, EAS3DVariantImpl> easVariant_;
    mutable Eigen::MatrixXd L;
    template <int enhancedStrainSize>
    void calculateDAndLMatrix(const auto& easFunction, const auto& par,
                              Eigen::Matrix<double, enhancedStrainSize, enhancedStrainSize>& DMat,
                              Eigen::MatrixXd& LMat) const {
      using namespace Dune;
      using namespace Dune::DerivativeDirections;

      auto strainFunction = DisplacementBasedElement::getStrainFunction(par);
      const auto C        = DisplacementBasedElement::getMaterialTangentFunction(par);
      const auto geo      = localView().element().geometry();
      const auto numNodes = DisplacementBasedElement::numberOfNodes;
      DMat.setZero();
      LMat.setZero(enhancedStrainSize, localView().size());
      for (const auto& [gpIndex, gp] : strainFunction.viewOverIntegrationPoints()) {
        const auto M      = easFunction.calcM(gp.position());
        const auto CEval  = C(gpIndex);
        const double detJ = geo.integrationElement(gp.position());
        DMat += M.transpose() * CEval * M * detJ * gp.weight();
        for (size_t i = 0U; i < numNodes; ++i) {
          const size_t I = Traits::worlddim * i;
          const auto Bi  = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(i)), on(gridElement));
          LMat.template block<enhancedStrainSize, Traits::worlddim>(0, I)
              += M.transpose() * CEval * Bi * detJ * gp.weight();
        }
      }
    }
  };
}  // namespace Ikarus

#else
#  error EnhancedAssumedStrains depends on dune-localfefunctions, which is not included
#endif
