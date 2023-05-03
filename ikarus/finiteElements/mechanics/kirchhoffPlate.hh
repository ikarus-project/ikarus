// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/localfefunctions/expressions/linearStrainsExpr.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feBases/scalarFE.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>

namespace Ikarus {

  template <typename Basis>
  class KirchhoffPlate : public Ikarus::ScalarFieldFE<typename Basis::FlatBasis> {
  public:
    using FlatBasis              = typename Basis::FlatBasis;
    using BaseDisp               = Ikarus::ScalarFieldFE<FlatBasis>;  // Handles globalIndices function
    using FERequirementType      = FErequirements<Eigen::VectorXd>;
    using ResultRequirementsType = ResultRequirements<Eigen::VectorXd>;
    using LocalView              = typename FlatBasis::LocalView;
    using GridView               = typename FlatBasis::GridView;

    KirchhoffPlate(const Basis& basis, const typename LocalView::Element& element, double p_Emodul, double p_nu,
                   double p_thickness)
        : BaseDisp(basis.flat(), element), Emodul{p_Emodul}, nu{p_nu}, thickness{p_thickness} {
      this->localView().bind(element);
      const int order = 2 * (this->localView().tree().finiteElement().localBasis().order());
      geometry_.emplace(this->localView().element().geometry());
      localBasis_ = Dune::CachedLocalBasis(this->localView().tree().finiteElement().localBasis());
      localBasis_.bind(Dune::QuadratureRules<double, Traits::mydim>::rule(this->localView().element().type(), order),
                       Dune::bindDerivatives(0, 1));
    }

    static Eigen::Matrix<double, 3, 3> constitutiveMatrix(double Emod, double p_nu, double p_thickness) {
      const double factor = Emod * Dune::power(p_thickness, 3) / (12.0 * (1.0 - p_nu * p_nu));
      Eigen::Matrix<double, 3, 3> D;
      D.setZero();
      D(0, 0) = 1;
      D(0, 1) = D(1, 0) = p_nu;
      D(1, 1)           = 1;
      D(2, 2)           = (1 - p_nu) / 2.0;
      D *= factor;
      return D;
    }

    using Traits = TraitsFromLocalView<LocalView>;

  public:
    double calculateScalar(const FERequirementType& par) const {
      const auto& wGlobal = par.getGlobalSolution(Ikarus::FESolutions::displacement);
      const auto& lambda  = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto D        = constitutiveMatrix(Emodul, nu, thickness);
      double energy       = 0.0;
      auto& ele           = this->localView().element();
      auto& fe            = this->localView().tree().finiteElement();
      Eigen::VectorX<double> wNodal;
      wNodal.setZero(fe.size());
      for (auto i = 0U; i < fe.size(); ++i)
        wNodal(i) = wGlobal[this->localView().index(this->localView().tree().localIndex(i))[0]];

      const auto& localBasis = fe.localBasis();

      const auto& rule = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * localBasis.order());
      /// Calculate Kirchhoff plate energy
      for (auto& gp : rule) {
        std::vector<Dune::FieldVector<double, 1>> dN_xixi;
        std::vector<Dune::FieldVector<double, 1>> dN_xieta;
        std::vector<Dune::FieldVector<double, 1>> dN_etaeta;
        std::vector<Dune::FieldVector<double, 1>> N_dune;
        Eigen::VectorXd N(fe.size());

        localBasis.evaluateFunction(gp.position(), N_dune);
        std::ranges::copy(N_dune, N.begin());
        localBasis.partial({2, 0}, gp.position(), dN_xixi);
        localBasis.partial({1, 1}, gp.position(), dN_xieta);
        localBasis.partial({0, 2}, gp.position(), dN_etaeta);

        const auto Jinv = toEigen(geometry_.jacobianInverseTransposed(gp.position())).transpose().eval();

        Eigen::VectorXd dN_xx(fe.size());
        Eigen::VectorXd dN_yy(fe.size());
        Eigen::VectorXd dN_xy(fe.size());
        using Dune::power;
        for (auto i = 0U; i < fe.size(); ++i) {
          dN_xx[i] = dN_xixi[i] * power(Jinv(0, 0), 2);
          dN_yy[i] = dN_etaeta[i] * power(Jinv(1, 1), 2);
          dN_xy[i] = dN_xieta[i] * Jinv(0, 0) * Jinv(1, 1);
        }
        Eigen::Vector<double, 3> kappa;
        kappa(0) = dN_xx.dot(wNodal);
        kappa(1) = dN_yy.dot(wNodal);
        kappa(2) = 2 * dN_xy.dot(wNodal);
        double w = N.dot(wNodal);

        energy += (0.5 * kappa.dot(D * kappa) - w * lambda) * geometry_.integrationElement(gp.position()) * gp.weight();
      }

      /// Clamp boundary using penalty method
      const double penaltyFactor = 1e8;
      if (ele.hasBoundaryIntersections())
        for (auto& intersection : intersections(this->localView().globalBasis().gridView(), ele))
          if (intersection.boundary()) {
            const auto& rule1 = Dune::QuadratureRules<double, 1>::rule(intersection.type(), 2 * localBasis.order());
            for (auto& gp : rule1) {
              const auto& gpInElement = intersection.geometryInInside().global(gp.position());
              std::vector<Dune::FieldMatrix<double, 1, 2>> dN_xi_eta;
              localBasis.evaluateJacobian(gpInElement, dN_xi_eta);
              Eigen::VectorXd dN_x(fe.size());
              Eigen::VectorXd dN_y(fe.size());
              const auto Jinv = toEigen(geometry_.jacobianInverseTransposed(gpInElement)).transpose().eval();
              for (auto i = 0U; i < fe.size(); ++i) {
                dN_x[i] = dN_xi_eta[i][0][0] * Jinv(0, 0);
                dN_y[i] = dN_xi_eta[i][0][1] * Jinv(1, 1);
              }
              const double w_x = dN_x.dot(wNodal);
              const double w_y = dN_y.dot(wNodal);

              energy += 0.0 * 0.5 * penaltyFactor * (w_x * w_x + w_y * w_y);
            }
          }

      return energy;
    }

    void calculateVector(const FERequirementType& par, typename Traits::VectorType& g) const {
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto D       = constitutiveMatrix(Emodul, nu, thickness);
      using namespace Dune::DerivativeDirections;
      auto& ele = this->localView().element();
      auto& fe  = this->localView().tree().finiteElement();

      const auto& localBasis = fe.localBasis();
      const auto geo         = this->localView().element().geometry();
      const auto& rule       = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * localBasis.order());

      g.template setZero(this->localView().size());
      for (const auto& gp : rule) {
        const auto Jinv         = toEigen(geo.jacobianInverseTransposed(gp.position())).transpose().eval();
        const double intElement = geo.integrationElement(gp.position()) * gp.weight();
        std::vector<Dune::FieldVector<double, 1>> shapeFunctionValues;
        localBasis.evaluateFunction(gp.position(), shapeFunctionValues);
        Eigen::VectorXd N;
        N.setZero(this->localView().size());
        for (size_t nn = 0; nn < this->localView().size(); ++nn)
          N[nn] = shapeFunctionValues[nn];
        for (size_t kk = 0; kk < this->localView().size(); ++kk)
          g[kk] -= N[kk] * lambda * intElement;
      }
    }

    void getShapeFunctionsAndDerivatives(const auto& pos,
                                         std::vector<Dune::FieldVector<double, 1>>& shapeFunctionValues,
                                         Eigen::Matrix2Xd& dNdX, Eigen::Matrix3Xd& ddNddX) const {
      auto& ele              = this->localView().element();
      auto& fe               = this->localView().tree().finiteElement();
      const auto geo         = this->localView().element().geometry();
      const auto& localBasis = fe.localBasis();
      const auto JInvT       = geo.jacobianInverseTransposed(pos);
      const auto ddXddXi     = toEigen(geo.impl().secondDerivativeOfPosition(pos)).eval();

      std::vector<Dune::FieldMatrix<double, 1, 2>> referenceGradients;
      localBasis.evaluateJacobian(pos, referenceGradients);
      std::vector<Dune::FieldVector<double, 2>> gradients(referenceGradients.size());

      for (size_t i = 0; i < gradients.size(); i++)
        JInvT.mv(referenceGradients[i][0], gradients[i]);
      localBasis.evaluateFunction(pos, shapeFunctionValues);

      dNdX.setZero(Eigen::NoChange, fe.size());
      ddNddX.setZero(Eigen::NoChange, fe.size());

      for (size_t i = 0; i < fe.size(); i++) {
        dNdX(0, i) = gradients[i][0];
        dNdX(1, i) = gradients[i][1];
      }

      using Dune::power;
      const auto Jacobian = toEigen(geo.jacobianTransposed(pos)).eval();
      Eigen::Matrix3d JJ;
      JJ << power(Jacobian(0, 0), 2), power(Jacobian(0, 1), 2), 2 * Jacobian(0, 0) * Jacobian(0, 1),
          power(Jacobian(1, 0), 2), power(Jacobian(1, 1), 2), 2 * Jacobian(1, 0) * Jacobian(1, 1),
          Jacobian(0, 0) * Jacobian(1, 0), Jacobian(0, 1) * Jacobian(1, 1),
          Jacobian(0, 0) * Jacobian(1, 1) + Jacobian(0, 1) * Jacobian(1, 0);
      auto JJInv = JJ.inverse();
      Eigen::Matrix3Xd ddNddXi;
      ddNddXi.setZero(Eigen::NoChange, fe.size());
      std::vector<Dune::FieldVector<double, 1>> dNdXiXi;
      std::vector<Dune::FieldVector<double, 1>> dNdEtaEta;
      std::vector<Dune::FieldVector<double, 1>> dNdXiEta;

      localBasis.partial({2, 0}, pos, dNdXiXi);
      localBasis.partial({0, 2}, pos, dNdEtaEta);
      localBasis.partial({1, 1}, pos, dNdXiEta);

      for (auto i = 0U; i < fe.size(); ++i) {
        ddNddXi(0, i) = dNdXiXi[i];
        ddNddXi(1, i) = dNdEtaEta[i];
        ddNddXi(2, i) = dNdXiEta[i];
      }
      ddNddX = JJInv * (ddNddXi - ddXddXi * dNdX);
      ddNddX.row(2) *= 2.0;
    }

    void calculateMatrix(const FERequirementType& par, typename Traits::MatrixType& h) const {
      using namespace Dune::DerivativeDirections;
      const auto& lambda = par.getParameter(Ikarus::FEParameter::loadfactor);
      const auto D       = constitutiveMatrix(Emodul, nu, thickness);
      auto& fe           = this->localView().tree().finiteElement();
      const auto& rule
          = Dune::QuadratureRules<double, 2>::rule(this->localView().element().type(), 2 * fe.localBasis().order());
      h.template setZero(this->localView().size(), this->localView().size());

      for (const auto& gp : rule) {
        const double intElement
            = this->localView().element().geometry().integrationElement(gp.position()) * gp.weight();
        std::vector<Dune::FieldVector<double, 1>> shapeFunctionValues;
        Eigen::Matrix2Xd dNdX;
        Eigen::Matrix3Xd ddNddX;
        getShapeFunctionsAndDerivatives(gp.position(), shapeFunctionValues, dNdX, ddNddX);
        Eigen::MatrixXd bop;
        bop.setZero(3, fe.size());
        for (auto i = 0U; i < fe.size(); ++i) {
          bop(0, i) = ddNddX(0, i);
          bop(1, i) = ddNddX(1, i);
          bop(2, i) = ddNddX(2, i);
        }
        h += bop.transpose() * D * bop * intElement;
      }
    }

    void calculateAt(const ResultRequirementsType& req, const Eigen::Vector<double, Traits::mydim>& local,
                     ResultTypeMap<double>& result) const {
      using namespace Dune::Indices;
      const auto& disp       = req.getGlobalSolution(Ikarus::FESolutions::displacement);
      const auto D           = constitutiveMatrix(Emodul, nu, thickness);
      auto& fe               = this->localView().tree().finiteElement();
      const auto& localBasis = fe.localBasis();
      const auto geo         = this->localView().element().geometry();
      auto gp                = toDune(local);
      Eigen::VectorXd local_disp;
      local_disp.setZero(this->localView().size());

      int disp_counter = 0;
      for (size_t i = 0; i < fe.size(); ++i) {
        auto globalIndex         = this->localView().index(this->localView().tree().localIndex(i));
        local_disp[disp_counter] = disp[globalIndex];
        disp_counter++;
      }

      std::vector<Dune::FieldVector<double, 1>> shapeFunctionValues;
      Eigen::Matrix2Xd dNdX;
      Eigen::Matrix3Xd ddNddX;
      getShapeFunctionsAndDerivatives(gp, shapeFunctionValues, dNdX, ddNddX);
      Eigen::MatrixXd bop;
      bop.setZero(3, fe.size());
      for (auto i = 0U; i < fe.size(); ++i) {
        bop(0, i) = ddNddX(0, i);
        bop(1, i) = ddNddX(1, i);
        bop(2, i) = ddNddX(2, i);
      }

      Eigen::Vector<double, 3> req_res;
      req_res.setZero();
      req_res = D * bop * local_disp;

      typename ResultTypeMap<double>::ResultArray resv;
      if (req.isResultRequested(ResultType::stressResultant)) {
        resv.setZero(3, 1);
        resv(0, 0) = req_res(0);
        resv(1, 0) = req_res(1);
        resv(2, 0) = req_res(2);
        result.insertOrAssignResult(ResultType::stressResultant, resv);
      }
    }

    Dune::CachedLocalBasis<std::remove_cvref_t<decltype(std::declval<LocalView>().tree().finiteElement().localBasis())>>
        localBasis_;

    /// Dune::Geometry<...> is not copy assignable, see https://gitlab.dune-project.org/core/dune-grid/-/issues/140,
    /// Thus, we wrap it inside a std::optional
    std::optional<typename LocalView::Element::Geometry> geometry_;

    double Emodul;
    double nu;
    double thickness;
  };

}  // namespace Ikarus