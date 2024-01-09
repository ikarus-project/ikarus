// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file kirchhoffloveshell.hh
 * @brief Definition of the KirchhoffLoveShell class for Kirchhoff-Love shell elements in Ikarus.
 */

#pragma once

#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <ikarus/finiteelements/febases/powerbasisfe.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/membranestrains.hh>
#include <ikarus/finiteelements/physicshelper.hh>

namespace Ikarus {

  /**
   * @brief Kirchhoff-Love shell finite element class.
   *
   * This class represents Kirchhoff-Love shell finite elements.
   *
   * @ingroup mechanics
   *
   * @tparam Basis_ The basis type for the finite element.
   * @tparam FERequirements_ The type representing the requirements for finite element calculations.
   * @tparam useEigenRef A boolean indicating whether to use Eigen references for efficiency.
   */
  template <typename Basis_, typename FERequirements_ = FERequirements<>, bool useEigenRef = false>
  class KirchhoffLoveShell : public PowerBasisFE<typename Basis_::FlatBasis> {
  public:
    using Basis                             = Basis_;
    using FlatBasis                         = typename Basis::FlatBasis;
    using BasePowerFE                       = PowerBasisFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType                 = FERequirements_;
    using ResultRequirementsType            = ResultRequirements<FERequirementType>;
    using LocalView                         = typename FlatBasis::LocalView;
    using Element                           = typename LocalView::Element;
    using Geometry                          = typename Element::Geometry;
    using GridView                          = typename FlatBasis::GridView;
    using Traits                            = TraitsFromLocalView<LocalView, useEigenRef>;
    static constexpr int myDim              = Traits::mydim;
    static constexpr int worlddim           = Traits::worlddim;
    static constexpr int membraneStrainSize = 3;
    static constexpr int bendingStrainSize  = 3;

    /**
     * \brief A structure representing kinematic variables.
     *
     * This structure holds various kinematic variables used in a mechanical analysis. It includes material tangent,
     * membrane strain, bending strain, Jacobian matrices of deformed and reference geometries, Hessian matrices of
     * deformed and reference geometries, the normal vector, and the normalized normal vector.
     *
     * \tparam ScalarType The scalar type for the matrix and vector elements.
     */
    template <typename ScalarType = double>
    struct KinematicVariables {
      Eigen::Matrix<double, 3, 3> C;      ///< material tangent
      Eigen::Vector3<ScalarType> epsV;    ///< membrane strain in Voigt notation
      Eigen::Vector3<ScalarType> kappaV;  ///< bending strain in Voigt notation
      Eigen::Matrix<ScalarType, 2, 3> j;  ///< Jacobian of the deformed geometry
      Eigen::Matrix<double, 2, 3> J;      ///< Jacobian of the reference geometry
      Eigen::Matrix3<ScalarType> h;       ///< Hessian of the deformed geometry
      Eigen::Matrix3<double> H;           ///< Hessian of the reference geometry
      Eigen::Vector3<ScalarType> a3N;     ///< Normal vector of the deformed geometry
      Eigen::Vector3<ScalarType> a3;      ///< normalized normal vector of the deformed geometry
    };

    /**
     * @brief Constructor for the KirchhoffLoveShell class.
     *
     * Initializes the KirchhoffLoveShell instance with the given parameters.
     *
     * @tparam VolumeLoad The type representing the volume load function.
     * @tparam NeumannBoundaryLoad The type representing the Neumann boundary load function.
     * @param globalBasis The global basis for the finite element.
     * @param element The local element to bind.
     * @param emod Young's modulus of the material.
     * @param nu Poisson's ratio of the material.
     * @param thickness Thickness of the shell.
     * @param p_volumeLoad The volume load function (optional, default is utils::LoadDefault).
     * @param p_neumannBoundary The Neumann boundary patch (optional, default is nullptr).
     * @param p_neumannBoundaryLoad The Neumann boundary load function (optional, default is LoadDefault).
     */
    template <typename VolumeLoad = utils::LoadDefault, typename NeumannBoundaryLoad = utils::LoadDefault>
    KirchhoffLoveShell(const Basis& globalBasis, const typename LocalView::Element& element, double emod, double nu,
                       double thickness, VolumeLoad p_volumeLoad = {},
                       const BoundaryPatch<GridView>* p_neumannBoundary = nullptr,
                       NeumannBoundaryLoad p_neumannBoundaryLoad        = {})
        : BasePowerFE(globalBasis.flat(), element),
          neumannBoundary{p_neumannBoundary},
          emod_{emod},
          nu_{nu},
          thickness_{thickness} {
      this->localView().bind(element);
      auto& first_child = this->localView().tree().child(0);
      const auto& fe    = first_child.finiteElement();
      numberOfNodes     = fe.size();
      dispAtNodes.resize(fe.size());
      order      = 2 * (fe.localBasis().order());
      localBasis = Dune::CachedLocalBasis(fe.localBasis());
      if constexpr (requires { this->localView().element().impl().getQuadratureRule(order); })
        if (this->localView().element().impl().isTrimmed())
          localBasis.bind(this->localView().element().impl().getQuadratureRule(order), Dune::bindDerivatives(0, 1, 2));
        else
          localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                          Dune::bindDerivatives(0, 1, 2));
      else
        localBasis.bind(Dune::QuadratureRules<double, myDim>::rule(this->localView().element().type(), order),
                        Dune::bindDerivatives(0, 1, 2));

      if constexpr (!std::is_same_v<VolumeLoad, utils::LoadDefault>) volumeLoad = p_volumeLoad;
      if constexpr (!std::is_same_v<NeumannBoundaryLoad, utils::LoadDefault>)
        neumannBoundaryLoad = p_neumannBoundaryLoad;

      assert(((not p_neumannBoundary and not neumannBoundaryLoad) or (p_neumannBoundary and neumannBoundaryLoad))
             && "If you pass a Neumann boundary you should also pass the function for the Neumann load!");
    }

  public:
    /**
     * @brief Get the displacement function and nodal displacements.
     *
     * Retrieves the displacement function and nodal displacements based on the given FERequirements.
     *
     * @tparam ScalarType The scalar type used for calculations.
     * @param par The FERequirements.
     * @param dx Optional additional displacement vector.
     * @return The displacement function.
     */
    template <typename ScalarType>
    auto getDisplacementFunction(const FERequirementType& par,
                                 const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      const auto& d = par.getGlobalSolution(FESolutions::displacement);
      auto geo      = std::make_shared<const typename GridView::GridView::template Codim<0>::Entity::Geometry>(
          this->localView().element().geometry());
      Dune::BlockVector<Dune::RealTuple<ScalarType, Traits::worlddim>> disp(dispAtNodes.size());

      if (dx)
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            disp[i][k2] = dx.value()[i * worlddim + k2]
                          + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
      else
        for (auto i = 0U; i < disp.size(); ++i)
          for (auto k2 = 0U; k2 < worlddim; ++k2)
            disp[i][k2] = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];

      Dune::StandardLocalFunction uFunction(localBasis, disp, geo);

      return uFunction;
    }

    /**
     * @brief Calculate the scalar value.
     *
     * Calculates the scalar value based on the given FERequirements.
     *
     * @param req The FERequirements.
     * @return The calculated scalar value.
     */
    double calculateScalar(const FERequirementType& req) const { return calculateScalarImpl<double>(req); }

    /**
     * @brief Calculate the vector associated with the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the calculation.
     * @param req The FERequirementType object specifying the requirements for the calculation.
     * @param force The vector to store the calculated result.
     */
    void calculateVector(const FERequirementType& req, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(req, force);
    }
    /**
     * @brief Calculate the matrix associated with the given FERequirementType.
     *
     * @tparam ScalarType The scalar type for the calculation.
     * @param req The FERequirementType object specifying the requirements for the calculation.
     * @param K The matrix to store the calculated result.
     */
    void calculateMatrix(const FERequirementType& req, typename Traits::template MatrixType<> K) const {
      calculateMatrixImpl<double>(req, K);
    }

    /**
     * @brief Calculate results at local coordinates.
     *
     * Calculates the results at the specified local coordinates based on the given requirements and stores them in the
     * result container.
     *
     * @param req The result requirements.
     * @param local The local coordinates at which results are to be calculated.
     * @param result The result container to store the calculated values.
     */
    void calculateAt([[maybe_unused]] const ResultRequirementsType& req,
                     [[maybe_unused]] const Dune::FieldVector<double, Traits::mydim>& local,
                     [[maybe_unused]] ResultTypeMap<double>& result) const {
      DUNE_THROW(Dune::NotImplemented, "No results are implemented");
    }

    Dune::CachedLocalBasis<
        std::remove_cvref_t<decltype(std::declval<LocalView>().tree().child(0).finiteElement().localBasis())>>
        localBasis;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Dune::FieldVector<double, Traits::worlddim>&,
                                                          const double&)>
        volumeLoad;
    std::function<Eigen::Vector<double, Traits::worlddim>(const Dune::FieldVector<double, Traits::worlddim>&,
                                                          const double&)>
        neumannBoundaryLoad;
    const BoundaryPatch<GridView>* neumannBoundary;
    mutable Dune::BlockVector<Dune::RealTuple<double, Traits::dimension>> dispAtNodes;
    DefaultMembraneStrain membraneStrain;
    double emod_;
    double nu_;
    double thickness_;
    size_t numberOfNodes{0};
    int order{};

  protected:
    /**
     * \brief Compute material properties and strains at a given integration point.
     *
     * \param gpPos The position of the integration point.
     * \param gpIndex The index of the integration point.
     * \param geo The geometry object providing position and derivatives.
     * \param uFunction The function representing the displacement field.
     *
     * \tparam gpPos The type of the integration point position.
     * \tparam gpIndex The type of the integration point index.
     * \tparam geo The type of the geometry object.
     * \tparam uFunction The type of the displacement field function.
     *
     * \return A tuple containing the material tangent, membrane strain, bending,
     *         Jacobian matrix of the reference position,  Jacobian matrix of the current position, Hessian matrix of
     * the current position, Hessian matrix of
     * the reference position, normal vector, and normalized normal vector at the given
     * integration point.
     */
    auto computeMaterialAndStrains(const Dune::FieldVector<double, 2>& gpPos, int gpIndex, const Geometry& geo,
                                   const auto& uFunction) const {
      using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

      KinematicVariables<ScalarType> kin;
      using namespace Dune;
      using namespace Dune::DerivativeDirections;
      const auto [X, Jd, Hd]              = geo.impl().zeroFirstAndSecondDerivativeOfPosition(gpPos);
      kin.J                               = toEigen(Jd);
      kin.H                               = toEigen(Hd);
      const Eigen::Matrix<double, 2, 2> A = kin.J * kin.J.transpose();
      Eigen::Matrix<double, 3, 3> G       = Eigen::Matrix<double, 3, 3>::Zero();

      G.block<2, 2>(0, 0)                    = A;
      G(2, 2)                                = 1;
      const Eigen::Matrix<double, 3, 3> GInv = G.inverse();
      kin.C                                  = materialTangent(GInv);

      kin.epsV = membraneStrain.value(gpPos, geo, uFunction);

      const auto& Ndd                             = localBasis.evaluateSecondDerivatives(gpIndex);
      const auto uasMatrix                        = Dune::viewAsEigenMatrixAsDynFixed(uFunction.coefficientsRef());
      const auto hessianu                         = Ndd.transpose().template cast<ScalarType>() * uasMatrix;
      kin.h                                       = kin.H + hessianu;
      const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(uFunction.evaluateDerivative(
          gpIndex, Dune::wrt(spatialAll), Dune::on(Dune::DerivativeDirections::referenceElement)));
      kin.j                                       = kin.J + gradu.transpose();
      kin.a3N                                     = (kin.j.row(0).cross(kin.j.row(1)));
      kin.a3                                      = kin.a3N.normalized();
      Eigen::Vector<ScalarType, 3> bV             = kin.h * kin.a3;
      bV(2) *= 2;  // Voigt notation requires the two here
      const auto BV = toVoigt(toEigen(geo.impl().secondFundamentalForm(gpPos)));
      kin.kappaV    = BV - bV;
      return kin;
    }

    template <typename ScalarType>
    void calculateMatrixImpl(const FERequirementType& par, typename Traits::template MatrixType<ScalarType> K,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      K.setZero();
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = getDisplacementFunction(par, dx);
      const auto& lambda   = par.getParameter(FEParameter::loadfactor);
      const auto geo       = this->localView().element().geometry();

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto intElement = geo.integrationElement(gp.position()) * gp.weight();
        const auto [C, epsV, kappaV, jE, J, h, H, a3N, a3]
            = computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);
        const Eigen::Vector<ScalarType, membraneStrainSize> membraneForces = thickness_ * C * epsV;
        const Eigen::Vector<ScalarType, bendingStrainSize> moments = Dune::power(thickness_, 3) / 12.0 * C * kappaV;

        const auto& Nd  = localBasis.evaluateJacobian(gpIndex);
        const auto& Ndd = localBasis.evaluateSecondDerivatives(gpIndex);
        for (size_t i = 0; i < numberOfNodes; ++i) {
          Eigen::Matrix<ScalarType, membraneStrainSize, worlddim> bopIMembrane
              = membraneStrain.derivative(gp.position(), jE, Nd, geo, uFunction, localBasis, i);

          Eigen::Matrix<ScalarType, bendingStrainSize, worlddim> bopIBending = bopBending(jE, h, Nd, Ndd, i, a3N, a3);
          for (size_t j = i; j < numberOfNodes; ++j) {
            auto KBlock = K.template block<worlddim, worlddim>(worlddim * i, worlddim * j);
            Eigen::Matrix<ScalarType, membraneStrainSize, worlddim> bopJMembrane
                = membraneStrain.derivative(gp.position(), jE, Nd, geo, uFunction, localBasis, j);
            Eigen::Matrix<ScalarType, bendingStrainSize, worlddim> bopJBending = bopBending(jE, h, Nd, Ndd, j, a3N, a3);
            KBlock += thickness_ * bopIMembrane.transpose() * C * bopJMembrane * intElement;
            KBlock += Dune::power(thickness_, 3) / 12.0 * bopIBending.transpose() * C * bopJBending * intElement;

            Eigen::Matrix<ScalarType, worlddim, worlddim> kgMembraneIJ
                = membraneStrain.secondDerivative(gp.position(), Nd, geo, uFunction, localBasis, membraneForces, i, j);
            Eigen::Matrix<ScalarType, worlddim, worlddim> kgBendingIJ
                = kgBending(jE, h, Nd, Ndd, a3N, a3, moments, i, j);
            KBlock += kgMembraneIJ * intElement;
            KBlock += kgBendingIJ * intElement;
          }
        }
      }
      K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
    }

    template <typename ScalarType>
    void calculateVectorImpl(const FERequirementType& par, typename Traits::template VectorType<ScalarType> force,
                             const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) const {
      force.setZero();
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = getDisplacementFunction(par, dx);
      const auto& lambda   = par.getParameter(FEParameter::loadfactor);
      const auto geo       = this->localView().element().geometry();

      // Internal forces
      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto [C, epsV, kappaV, jE, J, h, H, a3N, a3]
            = computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);
        const Eigen::Vector<ScalarType, 3> membraneForces = thickness_ * C * epsV;
        const Eigen::Vector<ScalarType, 3> moments        = Dune::power(thickness_, 3) / 12.0 * C * kappaV;

        const auto& Nd  = localBasis.evaluateJacobian(gpIndex);
        const auto& Ndd = localBasis.evaluateSecondDerivatives(gpIndex);
        for (size_t i = 0; i < numberOfNodes; ++i) {
          Eigen::Matrix<ScalarType, 3, 3> bopIMembrane
              = membraneStrain.derivative(gp.position(), jE, Nd, geo, uFunction, localBasis, i);
          Eigen::Matrix<ScalarType, 3, 3> bopIBending = bopBending(jE, h, Nd, Ndd, i, a3N, a3);
          force.template segment<3>(3 * i)
              += bopIMembrane.transpose() * membraneForces * geo.integrationElement(gp.position()) * gp.weight();
          force.template segment<3>(3 * i)
              += bopIBending.transpose() * moments * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // External forces volume forces over the domain
      if (volumeLoad) {
        const auto u = getDisplacementFunction(par, dx);
        for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
          Eigen::Vector<double, Traits::worlddim> fext = volumeLoad(geo.global(gp.position()), lambda);
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = uFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<worlddim>(worlddim * i)
                -= udCi * fext * geo.integrationElement(gp.position()) * gp.weight();
          }
        }
      }

      if (not neumannBoundary and not neumannBoundaryLoad) return;

      // External forces, boundary forces, i.e. at the Neumann boundary
      for (auto&& intersection : intersections(neumannBoundary->gridView(), this->localView().element())) {
        if (not neumannBoundary or not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = QuadratureRules<double, Element::mydimension - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const FieldVector<double, Element::mydimension>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double intElement = intersection.geometry().integrationElement(curQuad.position()) * curQuad.weight();
          for (size_t i = 0; i < numberOfNodes; ++i) {
            const auto udCi = uFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

            // Value of the Neumann data at the current position
            auto neumannValue = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda);
            force.template segment<worlddim>(worlddim * i) -= udCi * neumannValue * intElement;
          }
        }
      }
    }

    template <typename ScalarType>
    auto calculateScalarImpl(const FERequirementType& par, const std::optional<const Eigen::VectorX<ScalarType>>& dx
                                                           = std::nullopt) const -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto uFunction = getDisplacementFunction(par, dx);
      const auto& lambda   = par.getParameter(FEParameter::loadfactor);
      const auto geo       = this->localView().element().geometry();
      ScalarType energy    = 0.0;

      for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
        const auto [C, epsV, kappaV, j, J, h, H, a3N, a3]
            = computeMaterialAndStrains(gp.position(), gpIndex, geo, uFunction);

        const ScalarType membraneEnergy = 0.5 * thickness_ * epsV.dot(C * epsV);
        const ScalarType bendingEnergy  = 0.5 * Dune::power(thickness_, 3) / 12.0 * kappaV.dot(C * kappaV);
        energy += (membraneEnergy + bendingEnergy) * geo.integrationElement(gp.position()) * gp.weight();
      }

      if (volumeLoad) {
        for (const auto& [gpIndex, gp] : uFunction.viewOverIntegrationPoints()) {
          const auto u                                       = uFunction.evaluate(gpIndex);
          const Eigen::Vector<double, Traits::worlddim> fExt = volumeLoad(geo.global(gp.position()), lambda);
          energy -= u.dot(fExt) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      if (not neumannBoundary and not neumannBoundaryLoad) return energy;

      // line or surface loads, i.e., neumann boundary
      const auto& element = this->localView().element();
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary or not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = QuadratureRules<double, Traits::mydim - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const FieldVector<double, Traits::mydim>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double intElement = intersection.geometry().integrationElement(curQuad.position());

          // The value of the local function
          const auto u = uFunction.evaluate(quadPos);

          // Value of the Neumann data at the current position
          const auto neumannValue = neumannBoundaryLoad(intersection.geometry().global(curQuad.position()), lambda);
          energy -= neumannValue.dot(u) * curQuad.weight() * intElement;
        }
      }
      return energy;
    }

    template <typename ScalarType>
    Eigen::Matrix<ScalarType, 3, 3> kgBending(const Eigen::Matrix<ScalarType, 2, 3>& jcur,
                                              const Eigen::Matrix3<ScalarType>& h, const auto& dN, const auto& ddN,
                                              const Eigen::Vector3<ScalarType>& a3N,
                                              const Eigen::Vector3<ScalarType>& a3, const Eigen::Vector3<ScalarType>& S,
                                              int I, int J) const {
      Eigen::Matrix<ScalarType, 3, 3> kg;
      kg.setZero();

      const auto& dN1i = dN(I, 0);
      const auto& dN1j = dN(J, 0);
      const auto& dN2i = dN(I, 1);
      const auto& dN2j = dN(J, 1);

      const Eigen::Matrix<ScalarType, 3, 3> P
          = 1.0 / a3N.norm() * (Eigen::Matrix<double, 3, 3>::Identity() - a3 * a3.transpose());

      const auto a1dxI = Eigen::Matrix<double, 3, 3>::Identity()
                         * dN1i;  // the auto here allows the exploitation of the identity matrices,
      // due to Eigen's expression templates
      const auto a2dxI                            = Eigen::Matrix<double, 3, 3>::Identity() * dN2i;
      const auto a1dxJ                            = Eigen::Matrix<double, 3, 3>::Identity() * dN1j;
      const auto a2dxJ                            = Eigen::Matrix<double, 3, 3>::Identity() * dN2j;
      const auto a1                               = jcur.row(0);
      const auto a2                               = jcur.row(1);
      const Eigen::Matrix<ScalarType, 3, 3> a3NdI = a1dxI.colwise().cross(a2) - a2dxI.colwise().cross(a1);
      const Eigen::Matrix<ScalarType, 3, 3> a3NdJ = a1dxJ.colwise().cross(a2) - a2dxJ.colwise().cross(a1);
      Eigen::Matrix<ScalarType, 3, 3> a3dI        = P * a3NdI;
      Eigen::Matrix<ScalarType, 3, 3> a3dJ        = P * a3NdJ;
      for (int i = 0; i < 3; ++i) {
        const auto a_albe               = h.row(i).transpose();
        const auto& ddNI                = ddN(I, i);
        const auto& ddNJ                = ddN(J, i);
        Eigen::Vector3<ScalarType> vecd = P * a_albe;

        Eigen::Matrix<ScalarType, 3, 3> a3Ndd
            = 1.0 / a3N.squaredNorm()
              * ((3 * a3 * a3.transpose() - Eigen::Matrix<double, 3, 3>::Identity()) * (a3.dot(a_albe))
                 - a_albe * a3.transpose() - a3 * a_albe.transpose());

        Eigen::Matrix<ScalarType, 3, 3> secondDerivativeDirectorIJ = skew(((dN2i * dN1j - dN1i * dN2j) * vecd).eval());
        kg -= (a3NdI.transpose() * a3Ndd * a3NdJ + secondDerivativeDirectorIJ + (ddNI * a3dJ + ddNJ * a3dI.transpose()))
              * S[i] * (i == 2 ? 2 : 1);
      }

      return kg;
    }

    template <typename ScalarType>
    Eigen::Matrix<ScalarType, 3, 3> bopBending(const Eigen::Matrix<ScalarType, 2, 3>& jcur,
                                               const Eigen::Matrix3<ScalarType>& h, const auto& dN, const auto& ddN,
                                               const int node, const Eigen::Vector3<ScalarType>& a3N,
                                               const Eigen::Vector3<ScalarType>& a3) const {
      const Eigen::Matrix<ScalarType, 3, 3> a1dxI
          = Eigen::Matrix<double, 3, 3>::Identity() * dN(node, 0);  // this should be double
      // but the cross-product below complains otherwise
      const Eigen::Matrix<ScalarType, 3, 3> a2dxI = Eigen::Matrix<double, 3, 3>::Identity() * dN(node, 1);
      const auto a1                               = jcur.row(0);
      const auto a2                               = jcur.row(1);
      const Eigen::Matrix<ScalarType, 3, 3> a3NdI
          = a1dxI.colwise().cross(a2) - a2dxI.colwise().cross(a1);  // minus needed since order has
      // to be swapped to get column-wise cross product working
      const Eigen::Matrix<ScalarType, 3, 3> a3d1
          = 1.0 / a3N.norm() * (Eigen::Matrix<double, 3, 3>::Identity() - a3 * a3.transpose()) * a3NdI;

      Eigen::Matrix<ScalarType, 3, 3> bop = -(h * a3d1 + (a3 * ddN.row(node)).transpose());
      bop.row(2) *= 2;

      return bop;
    }

    Eigen::Matrix<double, 3, 3> materialTangent(const Eigen::Matrix<double, 3, 3>& Aconv) const {
      const double lambda   = emod_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
      const double mu       = emod_ / (2.0 * (1.0 + nu_));
      const double lambdbar = 2.0 * lambda * mu / (lambda + 2.0 * mu);
      Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>> moduli;
      const auto AconvT = tensorView(Aconv, std::array<Eigen::Index, 2>({3, 3}));
      moduli = lambdbar * dyadic(AconvT, AconvT).eval() + 2.0 * mu * symmetricFourthOrder<double>(Aconv, Aconv);

      auto C                          = toVoigt(moduli);
      Eigen::Matrix<double, 3, 3> C33 = C({0, 1, 5}, {0, 1, 5});
      return C33;
    }
  };
}  // namespace Ikarus
