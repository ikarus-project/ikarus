//
// Created by Alex on 26.06.2021.
//

#pragma once
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/utils/concepts.h>

#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"
#include <utility>

#include <dune/common/power.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Ikarus {

//  template <typename FEManager>
//  class ScalarAssembler {
//  public:
//    explicit ScalarAssembler(FEManager& dofManager) : feManager_{&dofManager} {}
//
//    double& getScalar(Ikarus::FiniteElements::ScalarAffordances scalarAffordances,
//                      const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      return getScalarImpl(scalarAffordances, feParameter);
//    }
//
//  private:
//    double& getScalarImpl(Ikarus::FiniteElements::ScalarAffordances scalarAffordances,
//                          const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      scal = 0.0;
//      for (auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
//        FiniteElements::FErequirements fErequirements;
//        fErequirements.scalarAffordances = scalarAffordances;
//        fErequirements.variables         = vars;
//        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
//        scal += calculateScalar(fe, fErequirements);
//      }
//      return scal;
//    }
//
//    FEManager* feManager_;
//    double scal{0.0};
//  };
//
//  template <typename FEManager, typename DirichletManager>
//  class VectorAssembler {
//  public:
//    explicit VectorAssembler(const FEManager& dofManager, const DirichletManager& dirichletManager)
//        : feManager_{&dofManager}, dirichletManager_{&dirichletManager} {}
//
//    Eigen::VectorXd& getVector(Ikarus::FiniteElements::VectorAffordances vectorAffordances,
//                               const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      return getVectorImpl(vectorAffordances, feParameter);
//    }
//
//    Eigen::VectorXd& getReducedVector(Ikarus::FiniteElements::VectorAffordances vectorAffordances,
//                                      const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      return getReducedVectorImpl(vectorAffordances, feParameter);
//    }
//
//    auto createFullVector(const Eigen::VectorXd& reducedVector) {
//      assert(reducedVector.size() == static_cast<Eigen::Index>(dirichletManager_->numberOfReducedDegreesOfFreedom())
//             && "The reduced vector you passed has the wrong dimensions.");
//      Eigen::Index reducedCounter = 0;
//      Eigen::VectorXd fullVec(feManager_->numberOfDegreesOfFreedom());
//      for (Eigen::Index i = 0; i < fullVec.size(); ++i) {
//        if (dirichletManager_->isConstrained(i)) {
//          ++reducedCounter;
//          fullVec[i] = 0.0;
//          continue;
//        } else
//          fullVec[i] = reducedVector[i - reducedCounter];
//      }
//      return fullVec;
//    }
//
//  private:
//    Eigen::VectorXd& getVectorImpl(Ikarus::FiniteElements::VectorAffordances vecAffordances,
//                                   const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      vec.setZero(feManager_->numberOfDegreesOfFreedom());
//      for (auto&& [fe, dofIndices, vars] : feManager_->elementIndicesVariableTuple()) {
//        FiniteElements::FErequirements fErequirements;
//        fErequirements.vectorAffordances = vecAffordances;
//        fErequirements.variables         = vars;
//        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
//        assert(dofIndices.size() == calculateVector(fe, fErequirements).size()
//               && "The returned vector has wrong rowSize!");
//        vec(dofIndices) += calculateVector(fe, fErequirements);
//      }
//      return vec;
//    }
//
//    Eigen::VectorXd& getReducedVectorImpl(Ikarus::FiniteElements::VectorAffordances vecAffordances,
//                                          const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      vecRed.setZero(dirichletManager_->numberOfReducedDegreesOfFreedom());
//      int reducedCounter = 0;
//      for (auto&& [fe, dofIndices, vars] : feManager_->elementIndicesVariableTuple()) {
//        FiniteElements::FErequirements fErequirements;
//        fErequirements.vectorAffordances = vecAffordances;
//        fErequirements.variables         = vars;
//        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
//
//        const auto f = calculateVector(fe, fErequirements);
//        assert(dofIndices.size() == f.size() && "The returned vector has wrong rowSize!");
//        for (int i = 0; auto&& dofIndex : dofIndices) {
//          if (dirichletManager_->isConstrained(dofIndex)) {
//            ++reducedCounter;
//            ++i;
//            continue;
//          } else
//            vecRed(dofIndex - dirichletManager_->constraintsBelow(dofIndex)) += f[i++];
//        }
//      }
//      return vecRed;
//    }
//
//    FEManager const* feManager_;
//    DirichletManager const* dirichletManager_;
//    Eigen::VectorXd vec{};
//    Eigen::VectorXd vecRed{};
//  };
//
//  template <typename FEManager, typename DirichletManager>
//  class DenseMatrixAssembler {
//  public:
//    explicit DenseMatrixAssembler(FEManager& dofManager, const DirichletManager& dirichletManager)
//        : feManager_{&dofManager}, dirichletManager_{&dirichletManager} {}
//
//    Eigen::MatrixXd& getMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances,
//                               const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      return getMatrixImpl(MatrixAffordances, feParameter);
//    }
//
//    Eigen::MatrixXd& getReducedMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances,
//                                      const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      return getReducedMatrixImpl(MatrixAffordances, feParameter);
//    }
//
//  private:
//    Eigen::MatrixXd& getMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matAffordances,
//                                   const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      mat.setZero(feManager_->numberOfDegreesOfFreedom(), feManager_->numberOfDegreesOfFreedom());
//      for (auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
//        FiniteElements::FErequirements fErequirements;
//        fErequirements.matrixAffordances = matAffordances;
//        fErequirements.variables         = vars;
//        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
//        assert(dofs.size() == calculateMatrix(fe, fErequirements).rows() && "The returned matrix has wrong rowSize!");
//        assert(dofs.size() == calculateMatrix(fe, fErequirements).cols() && "The returned matrix has wrong colSize!");
//        mat(dofs, dofs) += calculateMatrix(fe, fErequirements);
//      }
//      return mat;
//    }
//
//    Eigen::MatrixXd& getReducedMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matAffordances,
//                                          const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      matRed.setZero(dirichletManager_->numberOfReducedDegreesOfFreedom(),
//                     dirichletManager_->numberOfReducedDegreesOfFreedom());
//      for (auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
//        FiniteElements::FErequirements fErequirements;
//        fErequirements.matrixAffordances = matAffordances;
//        fErequirements.variables         = vars;
//        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
//        const auto eleMat = calculateMatrix(fe, fErequirements);
//        assert(dofs.size() == eleMat.rows() && "The returned matrix has wrong rowSize!");
//        assert(dofs.size() == eleMat.cols() && "The returned matrix has wrong colSize!");
//        std::cout<<"dofs: "<<dofs.transpose()<<std::endl;
//        for (int r = 0; r < dofs.size(); ++r) {
//          std::cout<<"r: "<<r<<" dofs[r]:"<<dofs[r];
//          if (dirichletManager_->isConstrained(dofs[r])) {
//            std::cout<<"skipped"<<std::endl;
//            continue;
//          } else {
//            for (int c = 0; c < dofs.size(); ++c) {
//              std::cout<<"c: "<<r<<" dofs[c]:"<<dofs[c];
//              if (dirichletManager_->isConstrained(dofs[c])) {
//                std::cout<<"skipped"<<std::endl;
//                continue;
//              }
//              std::cout<<"\ndirichletManager_->constraintsBelow(dofs[r]: "<<dirichletManager_->constraintsBelow(dofs[r]);
//              std::cout<<"\ndirichletManager_->constraintsBelow(dofs[c]: "<<dirichletManager_->constraintsBelow(dofs[c]);
//              std::cout<<"\nrrealindex"<<dofs[r] - dirichletManager_->constraintsBelow(dofs[r])<<std::endl;
//              std::cout<<"\ncrealindex"<<dofs[c] - dirichletManager_->constraintsBelow(dofs[c])<<std::endl;
//              matRed(dofs[r] - dirichletManager_->constraintsBelow(dofs[r]),
//                     dofs[c] - dirichletManager_->constraintsBelow(dofs[c]))
//                  += eleMat(r, c);
//            }
//          }
//        }
//      }
//      return matRed;
//    }
//
//    FEManager* feManager_;
//    DirichletManager const* dirichletManager_;
//    Eigen::MatrixXd mat{};
//    Eigen::MatrixXd matRed{};
//  };
//
//  template <typename FEManager, typename DirichletManager>
//  class SparseMatrixAssembler {
//  public:
//
//    explicit SparseMatrixAssembler(FEManager& dofManager, const DirichletManager& dirichletManager)
//        : feManager_{&dofManager}, dirichletManager_{&dirichletManager} {}
//
//    using GridView = typename FEManager::GridView;
//
//    Eigen::SparseMatrix<double>& getMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances,
//                                           const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      return getMatrixImpl(MatrixAffordances, feParameter);
//    }
//
//    Eigen::SparseMatrix<double>& getReducedMatrix(Ikarus::FiniteElements::MatrixAffordances MatrixAffordances,
//                                                  const std::optional<FEParameterValuePair>& feParameter
//                                                  = std::nullopt) {
//      return getReducedMatrixImpl(MatrixAffordances, feParameter);
//    }
//
//  private:
//    Eigen::SparseMatrix<double>& getMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matrixAffordances,
//                                               const std::optional<FEParameterValuePair>& feParameter = std::nullopt) {
//      if (!isOccupationPatternCreated) createOccupationPattern();
//      if (!arelinearDofsPerElementCreated) createlinearDofsPerElement();
//      spMat.coeffs().setZero();
//      Eigen::MatrixXd A;
//      for (size_t elementIndex = 0; auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
//        FiniteElements::FErequirements fErequirements;
//        fErequirements.matrixAffordances = matrixAffordances;
//        fErequirements.variables         = vars;
//        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
//        A = calculateMatrix(fe, fErequirements);
//        assert(dofs.size() == A.rows() && "The returned matrix has wrong rowSize!");
//        assert(dofs.size() == A.cols() && "The returned matrix has wrong colSize!");
//        for (Eigen::Index linearIndex = 0; double matrixEntry : A.reshaped())
//          spMat.coeffs()(elementLinearIndices[elementIndex][linearIndex++]) += matrixEntry;
//        ++elementIndex;
//      }
//      return spMat;
//    }
//
//    Eigen::SparseMatrix<double>& getReducedMatrixImpl(Ikarus::FiniteElements::MatrixAffordances matrixAffordances,
//                                                      const std::optional<FEParameterValuePair>& feParameter
//                                                      = std::nullopt) {
//      if (!isReducedOccupationPatternCreated) createReducedOccupationPattern();
//      if (!arelinearReducedDofsPerElementCreated) createlinearDofsPerElementReduced();
//      spMatReduced.coeffs().setZero();
//      Eigen::MatrixXd A;
//      for (size_t elementIndex = 0; auto [fe, dofs, vars] : feManager_->elementIndicesVariableTuple()) {
//        FiniteElements::FErequirements fErequirements;
//        fErequirements.matrixAffordances = matrixAffordances;
//        fErequirements.variables         = vars;
//        if (feParameter) fErequirements.parameter.insert({feParameter.value().type, feParameter.value().value});
//        A = calculateMatrix(fe, fErequirements);
//        assert(dofs.size() == A.rows() && "The returned matrix has wrong rowSize!");
//        assert(dofs.size() == A.cols() && "The returned matrix has wrong colSize!");
//        Eigen::Index linearIndex = 0;
//        for (int r = 0; r < dofs.size(); ++r) {
//          if (dirichletManager_->isConstrained(dofs[r]))
//            continue;
//          else {
//            for (int c = 0; c < dofs.size(); ++c) {
//              if (dirichletManager_->isConstrained(dofs[c])) continue;
//              spMatReduced.coeffs()(elementLinearReducedIndices[elementIndex][linearIndex++]) += A(r, c);
//            }
//          }
//        }
//        ++elementIndex;
//      }
//      return spMatReduced;
//    }
//
//    // https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
//    void createOccupationPattern() {
//      spMat.resize(feManager_->numberOfDegreesOfFreedom(), feManager_->numberOfDegreesOfFreedom());
//      std::vector<Eigen::Triplet<double>> vectorOfTriples;
//      int estimateOfConnectivity = 8;
//      using std::size;
//      vectorOfTriples.reserve(estimateOfConnectivity * feManager_->getGridView()->size(GridView::dimension));
//      for (auto&& dofsOfElement : feManager_->elementDofs())
//        for (auto&& c : dofsOfElement)
//          for (auto&& r : dofsOfElement)
//            vectorOfTriples.emplace_back(r, c, 0.0);
//
//      spMat.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
//      isOccupationPatternCreated = true;
//    }
//
//    // https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
//    void createReducedOccupationPattern() {
//      spMatReduced.resize(dirichletManager_->numberOfReducedDegreesOfFreedom(),
//                          dirichletManager_->numberOfReducedDegreesOfFreedom());
//      std::vector<Eigen::Triplet<double>> vectorOfTriples;
//      const int estimateOfConnectivity = 8;
//      using std::size;
//
//      vectorOfTriples.reserve(estimateOfConnectivity * feManager_->getGridView()->size(GridView::dimension));
//      for (auto&& dofs : feManager_->elementDofs()) {
//        for (int r = 0; r < dofs.size(); ++r) {
//          if (dirichletManager_->isConstrained(dofs[r]))
//            continue;
//          else {
//            for (int c = 0; c < dofs.size(); ++c) {
//              if (dirichletManager_->isConstrained(dofs[c])) continue;
//              vectorOfTriples.emplace_back(dofs[r] - dirichletManager_->constraintsBelow(dofs[r]),
//                                           dofs[c] - dirichletManager_->constraintsBelow(dofs[c]), 0.0);
//            }
//          }
//        }
//      }
//
//      spMatReduced.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
//      isReducedOccupationPatternCreated = true;
//    }
//
//    // This function save the indices of each element in the underlying vector which stores the sparse matrix entries
//    void createlinearDofsPerElement() {
//      for (Eigen::Index eleIndex = 0; auto&& dofsOfElement : feManager_->elementDofs()) {
//        elementLinearIndices.emplace_back(Dune::Power<2>::eval(dofsOfElement.size()));
//        for (Eigen::Index linearIndexOfElement = 0; auto&& c : dofsOfElement)
//          for (auto&& r : dofsOfElement)
//            elementLinearIndices.back()[linearIndexOfElement++] = spMat.getLinearIndex(r, c);
//        ++eleIndex;
//      }
//      arelinearDofsPerElementCreated = true;
//    }
//
//    // This function save the indices of each element in the underlying vector which stores the sparse matrix entries
//    void createlinearDofsPerElementReduced() {
//      for (auto&& dofs : feManager_->elementDofs()) {
//        elementLinearReducedIndices.emplace_back();
//        for (int r = 0; r < dofs.size(); ++r) {
//          if (dirichletManager_->isConstrained(dofs[r])) continue;
//          for (int c = 0; c < dofs.size(); ++c) {
//            if (dirichletManager_->isConstrained(dofs[c])) continue;
//            elementLinearReducedIndices.back().push_back(
//                spMatReduced.getLinearIndex(dofs[r] - dirichletManager_->constraintsBelow(dofs[r]),
//                                            dofs[c] - dirichletManager_->constraintsBelow(dofs[c])));
//          }
//        }
//      }
//      arelinearReducedDofsPerElementCreated = true;
//    }
//
//    Eigen::SparseMatrix<double> spMat;
//    Eigen::SparseMatrix<double> spMatReduced;
//    bool isOccupationPatternCreated{false};
//    bool isReducedOccupationPatternCreated{false};
//    bool arelinearDofsPerElementCreated{false};
//    bool arelinearReducedDofsPerElementCreated{false};
//    FEManager* feManager_;
//    DirichletManager const* dirichletManager_;
//    std::vector<std::vector<Eigen::Index>> elementLinearIndices;
//    std::vector<std::vector<Eigen::Index>> elementLinearReducedIndices;
//  };


  template <typename Basis, typename FEContainer> //requires Ikarus::Concepts::FlatIndexBasis<Basis>
  class DenseFlatAssembler {
  public:
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    explicit DenseFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {}

    Eigen::MatrixXd& getMatrix(const Ikarus::MatrixAffordances& p_matrixAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
      return getMatrixImpl(p_matrixAffordances,displacement, lambda);
    }

    Eigen::MatrixXd& getMatrix(const Ikarus::MatrixAffordances& p_matrixAffordances) {
      return getMatrixImpl(p_matrixAffordances);
    }

    Eigen::VectorXd& getVector(const Ikarus::VectorAffordances& p_vectorAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
      return getVectorImpl(p_vectorAffordances,displacement, lambda);
    }

    Eigen::VectorXd& getVector(const Ikarus::VectorAffordances& p_vectorAffordances) {
      return getVectorImpl(p_vectorAffordances);
    }

    double getScalar(const Ikarus::ScalarAffordances & p_scalarAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
      return getScalarImpl(p_scalarAffordances,displacement, lambda);
    }

    double getScalar(const Ikarus::ScalarAffordances & p_scalarAffordances) {
      return getScalarImpl(p_scalarAffordances);
    }

  private:
    Eigen::MatrixXd& getMatrixImpl(const Ikarus::MatrixAffordances& p_matrixAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
      mat.setZero(basis_->size(), basis_->size());
      RequirementType requirements;
      requirements.matrixAffordances = p_matrixAffordances;
      requirements.sols.emplace_back(displacement);
      requirements.parameter.insert({Ikarus::FEParameter::loadfactor,lambda});
      for (auto& fe : feContainer) {
        auto matLoc      = fe.calculateMatrix(requirements);
        auto globalIndices = fe.globalIndices();
        for (auto i = 0; auto idi : fe.globalIndices()) {
          for (auto j = 0; auto idj : fe.globalIndices()) {
            mat(idi[0], idj[0]) += matLoc(i, j);
            ++j;
          }
          ++i;
        }
      }
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) mat.col(i).setZero();
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) mat.row(i).setZero();
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) mat(i, i) = 1;
      return mat;
    }

    Eigen::MatrixXd& getMatrixImpl(const Ikarus::MatrixAffordances& p_matrixAffordances) {
      mat.setZero(basis_->size(), basis_->size());
      RequirementType requirements;
      requirements.matrixAffordances = p_matrixAffordances;
      for (auto& fe : feContainer) {
        auto matLoc      = fe.calculateMatrix(requirements);
        auto globalIndices = fe.globalIndices();
        for (auto i = 0; auto idi : fe.globalIndices()) {
          for (auto j = 0; auto idj : fe.globalIndices()) {
            mat(idi[0], idj[0]) += matLoc(i, j);
            ++j;
          }
          ++i;
        }
      }
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) mat.col(i).setZero();
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) mat.row(i).setZero();
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) mat(i, i) = 1;
      return mat;
    }

    Eigen::VectorXd& getVectorImpl(const Ikarus::VectorAffordances& p_vectorAffordances, const Eigen::VectorXd& displacement, const double& lambda) {
      vec.setZero(basis_->size());
      auto localView = basis_->localView();
      RequirementType requirements;
      requirements.vectorAffordances = p_vectorAffordances;
      requirements.sols.emplace_back(displacement);
      requirements.parameter.insert({Ikarus::FEParameter::loadfactor,lambda});
      for (auto& fe : feContainer) {
        //      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)> fe(localView, 1000, 0.3);

        auto vecLocal = fe.calculateVector(requirements);
        for (int i = 0;auto id : fe.globalIndices()) {
          vec(id[0]) += vecLocal(i);
          ++i;
        }
      }
      for (auto i = 0U; i < basis_->size(); ++i) {
        if (dirichletFlags->at(i)) vec[i] = 0;
      }

      return vec;
    }

    Eigen::VectorXd& getVectorImpl(const Ikarus::VectorAffordances& p_vectorAffordances) {
      vec.setZero(basis_->size());
      auto localView = basis_->localView();
      RequirementType requirements;
      requirements.vectorAffordances = p_vectorAffordances;
      for (auto& fe : feContainer) {
        auto vecLocal = fe.calculateVector(requirements);
        for (int i = 0;auto id : fe.globalIndices()) {
          vec(id[0]) += vecLocal(i);
          ++i;
        }
      }
      for (auto i = 0U; i < basis_->size(); ++i) {
        if (dirichletFlags->at(i)) vec[i] = 0;
      }

      return vec;
    }

    double getScalarImpl(const Ikarus::ScalarAffordances& p_scalarAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
      double scalar = 0.0;
      vec.setZero(basis_->size());
      RequirementType requirements;
      requirements.scalarAffordances = p_scalarAffordances;
      requirements.sols.emplace_back(displacement);
      requirements.parameter.insert({Ikarus::FEParameter::loadfactor,lambda});
      auto localView = basis_->localView();

      for (auto& fe : feContainer) {
        //      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)> fe(localView, 1000, 0.3);

        for (int i = 0;auto id : fe.globalIndices()) {
          scalar += fe.calculateScalar(requirements);
          ++i;
        }
      }

      return scalar;
    }

    double getScalarImpl(const Ikarus::ScalarAffordances& p_scalarAffordances) {
      double scalar = 0.0;
      vec.setZero(basis_->size());
      RequirementType requirements;
      requirements.scalarAffordances = p_scalarAffordances;
      auto localView = basis_->localView();

      for (auto& fe : feContainer) {
        for (int i = 0;auto id : fe.globalIndices()) {
          scalar += fe.calculateScalar(requirements);
          ++i;
        }
      }

      return scalar;
    }

    Basis const* basis_;
    FEContainer const& feContainer;
    std::vector<bool> const* dirichletFlags;
    Eigen::MatrixXd mat{};
    Eigen::VectorXd vec{};
  };

}  // namespace Ikarus::Assembler