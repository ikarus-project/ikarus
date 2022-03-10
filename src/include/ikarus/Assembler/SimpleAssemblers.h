//
// Created by Alex on 26.06.2021.
//

#pragma once
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/utils/concepts.h>

#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"
#include <ranges>
#include <utility>

#include <dune/common/power.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Ikarus {

//  namespace Impl {
//    template <typename FEContainer>
//    requires requires { {std::declval<typename FEContainer::value_type>().globalIndices(std::declval<std::vector<typename FEContainer::value_type::GlobalIndex>&>())}; }
//    auto dofsOfElements(const FEContainer& feContainer) {
//      std::vector<typename FEContainer::value_type::GlobalIndex> dofs;
//      return feContainer | std::views::transform([dofs = move(dofs)]  (auto&& fe)mutable {
//               dofs.resize(0);
//               fe.globalIndices(dofs);
//               return dofs; });
//    }
//  }  // namespace Impl
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
  //        assert(dofs.size() == calculateMatrix(fe, fErequirements).rows() && "The returned matrix has wrong
  //        rowSize!"); assert(dofs.size() == calculateMatrix(fe, fErequirements).cols() && "The returned matrix has
  //        wrong colSize!"); mat(dofs, dofs) += calculateMatrix(fe, fErequirements);
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
  //              std::cout<<"\ndirichletManager_->constraintsBelow(dofs[r]:
  //              "<<dirichletManager_->constraintsBelow(dofs[r]);
  //              std::cout<<"\ndirichletManager_->constraintsBelow(dofs[c]:
  //              "<<dirichletManager_->constraintsBelow(dofs[c]); std::cout<<"\nrrealindex"<<dofs[r] -
  //              dirichletManager_->constraintsBelow(dofs[r])<<std::endl; std::cout<<"\ncrealindex"<<dofs[c] -
  //              dirichletManager_->constraintsBelow(dofs[c])<<std::endl; matRed(dofs[r] -
  //              dirichletManager_->constraintsBelow(dofs[r]),
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

  template <typename Basis, typename FEContainer>  // requires Ikarus::Concepts::FlatIndexBasis<BasisEmbedded>
  class SparseFlatAssembler {
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    using GlobalIndex = typename FEContainer::value_type::GlobalIndex;

  public:
    SparseFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {
      constraintsBelow_.reserve(basis_->size());
      size_t counter = 0;
      for (auto iv : std::ranges::iota_view{size_t(0), basis_->size()}) {
        constraintsBelow_.emplace_back(counter);
        if (dirichFlags[iv]) ++counter;
      }

      fixedDofs = std::ranges::count(dirichFlags, true);
    }

    using GridView = typename Basis::GridView;

    Eigen::SparseMatrix<double>& getMatrix(const RequirementType& fErequirements) {
      return getMatrixImpl(fErequirements);
    }

    Eigen::SparseMatrix<double>& getReducedMatrix(const RequirementType& fErequirements) {
      return getReducedMatrixImpl(fErequirements);
    }

  private:
    Eigen::SparseMatrix<double>& getMatrixImpl(const RequirementType& fErequirements) {
      if (!isOccupationPatternCreated) createOccupationPattern();
      if (!arelinearDofsPerElementCreated) createlinearDofsPerElement();
      spMat.coeffs().setZero();
      Eigen::MatrixXd A;
      for (size_t elementIndex = 0; const auto& fe : feContainer) {
        A.setZero(fe.size(),fe.size());
        fe.calculateMatrix(fErequirements,A);
        assert(std::sqrt(elementLinearIndices[elementIndex].size()) == A.rows()
               && "The returned matrix has wrong rowSize!");
        assert(std::sqrt(elementLinearIndices[elementIndex].size()) == A.cols()
               && "The returned matrix has wrong colSize!");
        for (Eigen::Index linearIndex = 0; double matrixEntry : A.reshaped())
          spMat.coeffs()(elementLinearIndices[elementIndex][linearIndex++]) += matrixEntry;
        ++elementIndex;
      }
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) spMat.col(i) *= 0;
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) spMat.row(i) *= 0;
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) spMat.diagonal()[i] = 1;
      return spMat;
    }

    Eigen::SparseMatrix<double>& getReducedMatrixImpl(const RequirementType& fErequirements) {
      if (!isReducedOccupationPatternCreated) createReducedOccupationPattern();
      if (!arelinearReducedDofsPerElementCreated) createlinearDofsPerElementReduced();
      spMatReduced.coeffs().setZero();
      Eigen::MatrixXd A;
      std::vector<GlobalIndex> dofs;
      for (size_t elementIndex = 0; const auto& fe : feContainer) {
        A.setZero(fe.size(),fe.size());
        dofs.resize(0);
        fe.calculateMatrix(fErequirements,A);
        fe.globalIndices(dofs);
        assert(dofs.size() == static_cast<unsigned>(A.rows()) && "The returned matrix has wrong rowSize!");
        assert(dofs.size() == static_cast<unsigned>(A.cols()) && "The returned matrix has wrong colSize!");
        Eigen::Index linearIndex = 0;
        for (auto r = 0U; r < dofs.size(); ++r) {
          if (dirichletFlags->at(dofs[r][0]))
            continue;
          else {
            for (auto c = 0U; c < dofs.size(); ++c) {
              if (dirichletFlags->at(dofs[c][0])) continue;
              spMatReduced.coeffs()(elementLinearReducedIndices[elementIndex][linearIndex++]) += A(r, c);
            }
          }
        }
        ++elementIndex;
      }
      return spMatReduced;
    }

    // https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
    void createOccupationPattern() {
      spMat.resize(basis_->size(), basis_->size());
      std::vector<Eigen::Triplet<double>> vectorOfTriples;
      int estimateOfConnectivity = 8;
      using std::size;
      vectorOfTriples.reserve(estimateOfConnectivity * basis_->gridView().size(GridView::dimension));
      std::vector<GlobalIndex> dofs;
      for (auto&& fe : feContainer) {
        dofs.resize(0);
        fe.globalIndices(dofs);
        for (auto idi : dofs)
          for (auto idj : dofs)
            vectorOfTriples.emplace_back(idi[0], idj[0], 0.0);
      }

      spMat.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
      isOccupationPatternCreated = true;
    }

    // https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
    void createReducedOccupationPattern() {
      spMatReduced.resize(reducedSize(), reducedSize());
      std::vector<Eigen::Triplet<double>> vectorOfTriples;
      const int estimateOfConnectivity = 8;
      using std::size;

      vectorOfTriples.reserve(estimateOfConnectivity * basis_->gridView().size(GridView::dimension));
      std::vector<GlobalIndex> dofs;
      for (auto& fe : feContainer) {
        dofs.resize(0);
        fe.globalIndices(dofs);
        for (auto r = 0U; r < dofs.size(); ++r) {
          if (dirichletFlags->at(dofs[r][0]))
            continue;
          else {
            for (auto c = 0U; c < dofs.size(); ++c) {
              if (dirichletFlags->at(dofs[c][0])) continue;
              vectorOfTriples.emplace_back(dofs[r][0] - constraintsBelow_[dofs[r][0]], dofs[c][0] - constraintsBelow_[dofs[c][0]],
                                           0.0);
            }
          }
        }
      }

      spMatReduced.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
      isReducedOccupationPatternCreated = true;
    }

    // This function save the indices of each element in the underlying vector which stores the sparse matrix entries
    void createlinearDofsPerElement() {
      std::vector<GlobalIndex> dofs;
      for (auto&& fe : feContainer) {
        dofs.resize(0);
        fe.globalIndices(dofs);
        elementLinearIndices.emplace_back(Dune::Power<2>::eval(dofs.size()));
        for (Eigen::Index linearIndexOfElement = 0; auto&& c : dofs)
          for (auto&& r : dofs)
            elementLinearIndices.back()[linearIndexOfElement++] = spMat.getLinearIndex(r[0], c[0]);
      }
      arelinearDofsPerElementCreated = true;
    }

    // This function save the indices of each element in the underlying vector which stores the sparse matrix entries
    void createlinearDofsPerElementReduced() {
      std::vector<GlobalIndex> dofs;
      for (auto&& fe : feContainer) {
        dofs.resize(0);
        fe.globalIndices(dofs);
        elementLinearReducedIndices.emplace_back();
        for (auto r = 0U; r < dofs.size(); ++r) {
          if (dirichletFlags->at(dofs[r][0])) continue;
          for (auto c = 0U; c < dofs.size(); ++c) {
            if (dirichletFlags->at(dofs[c][0])) continue;
            elementLinearReducedIndices.back().push_back(spMatReduced.getLinearIndex(
                dofs[r][0] - constraintsBelow_[dofs[r][0]], dofs[c][0] - constraintsBelow_[dofs[c][0]]));
          }
        }
      }
      arelinearReducedDofsPerElementCreated = true;
    }

    size_t reducedSize(int i = 0) { return basis_->size() - fixedDofs; }

    size_t size(int i = 0) { return basis_->size(); }

  private:
    Eigen::SparseMatrix<double> spMat;
    Eigen::SparseMatrix<double> spMatReduced;
    bool isOccupationPatternCreated{false};
    bool isReducedOccupationPatternCreated{false};
    bool arelinearDofsPerElementCreated{false};
    bool arelinearReducedDofsPerElementCreated{false};
    Basis const* basis_;
    FEContainer const& feContainer;
    std::vector<bool> const* dirichletFlags;
    std::vector<size_t> constraintsBelow_;
    std::vector<std::vector<Eigen::Index>> elementLinearIndices;
    std::vector<std::vector<Eigen::Index>> elementLinearReducedIndices;
    size_t fixedDofs{};
  };

  template <typename Basis, typename FEContainer>
  class VectorFlatAssembler {
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    using GlobalIndex = typename FEContainer::value_type::GlobalIndex;
  public:
    VectorFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {
      constraintsBelow_.reserve(basis->size());
      size_t counter = 0;
      for (auto iv : std::ranges::iota_view{size_t(0), basis->size()}) {
        if (dirichFlags[iv]) ++counter;
        constraintsBelow_.emplace_back(counter--);
      }
      fixedDofs = std::ranges::count(dirichFlags, true);
    }

    Eigen::VectorXd& getVector(const RequirementType& fErequirements) { return getVectorImpl(fErequirements); }

    Eigen::VectorXd& getReducedVector(const RequirementType& fErequirements) {
      return getReducedVectorImpl(fErequirements);
    }

    auto createFullVector(const Eigen::VectorXd& reducedVector) {
      assert(reducedVector.size() == static_cast<Eigen::Index>(this->reducedSize())
             && "The reduced vector you passed has the wrong dimensions.");
      Eigen::Index reducedCounter = 0;
      Eigen::VectorXd fullVec(basis_->size());
      for (Eigen::Index i = 0; i < fullVec.size(); ++i) {
        if (dirichletFlags->at(i)) {
          ++reducedCounter;
          fullVec[i] = 0.0;
          continue;
        } else
          fullVec[i] = reducedVector[i - reducedCounter];
      }
      return fullVec;
    }

    size_t reducedSize(int i = 0) { return basis_->size() - fixedDofs; }

    size_t size(int i = 0) { return basis_->size(); }

  private:
    Eigen::VectorXd& getVectorImpl(const RequirementType& fErequirements) {
      vec.setZero(basis_->size());
      Eigen::VectorXd vecLocal;
      std::vector<GlobalIndex> dofs;
      for (auto& fe : feContainer) {
        dofs.resize(0);
        vecLocal.setZero(fe.size());
        fe.globalIndices(dofs);
        fe.calculateVector(fErequirements,vecLocal);
        vec(dofs) += vecLocal;
      }

      return vec;
    }

    Eigen::VectorXd& getReducedVectorImpl(const RequirementType& fErequirements) {
      vecRed.setZero(reducedSize());
      int reducedCounter = 0;
      for (auto& fe : feContainer) {
        const auto f    = fe.calculateVector(fErequirements);
        const auto dofs = f.globalIndices();
        assert(dofs.size() == f.size() && "The returned vector has wrong rowSize!");
        for (int i = 0; auto&& dofIndex : dofs) {
          if (dirichletFlags[dofIndex]) {
            ++reducedCounter;
            ++i;
            continue;
          } else
            vecRed(dofIndex[0] - constraintsBelow_[dofIndex[0]]) += f[i++];
        }
      }
      return vecRed;
    }

    Basis const* basis_;
    FEContainer const& feContainer;
    std::vector<bool> const* dirichletFlags;
    std::vector<size_t> constraintsBelow_;
    Eigen::VectorXd vec{};
    Eigen::VectorXd vecRed{};
    size_t fixedDofs{};
  };

  template <typename Basis, typename FEContainer>
  class ScalarAssembler {
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    friend SparseFlatAssembler<Basis, FEContainer>;
    friend VectorFlatAssembler<Basis, FEContainer>;

  public:
    ScalarAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {}

    double& getScalar(const RequirementType& fErequirements) { return getScalarImpl(fErequirements); }

  private:
    double& getScalarImpl(const RequirementType& fErequirements) {
      scal = 0.0;
      for (auto& fe : feContainer)
        scal += fe.calculateScalar(fErequirements);
      return scal;
    }

    Basis const* basis_;
    FEContainer const& feContainer;
    std::vector<bool> const* dirichletFlags;
    double scal{0.0};
  };

  template <typename Basis, typename FEContainer>  // requires Ikarus::Concepts::FlatIndexBasis<BasisEmbedded>
  class DenseFlatSimpleAssembler {
  public:
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    using GlobalIndex = typename FEContainer::value_type::GlobalIndex;
    explicit DenseFlatSimpleAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {}

    Eigen::MatrixXd& getMatrix(const Ikarus::MatrixAffordances& p_matrixAffordances,
                               const Eigen::VectorXd& displacement, const double& lambda) {
      return getMatrixImpl(p_matrixAffordances, displacement, lambda);
    }

    Eigen::VectorXd& getVector(const Ikarus::VectorAffordances& p_vectorAffordances,
                               const Eigen::VectorXd& displacement, const double& lambda) {
      return getVectorImpl(p_vectorAffordances, displacement, lambda);
    }

    double getScalar(const Ikarus::ScalarAffordances& p_scalarAffordances, const Eigen::VectorXd& displacement,
                     const double& lambda) {
      return getScalarImpl(p_scalarAffordances, displacement, lambda);
    }

  private:
    Eigen::MatrixXd& getMatrixImpl(const Ikarus::MatrixAffordances& p_matrixAffordances,
                                   const Eigen::VectorXd& displacement, const double& lambda) {
      mat.setZero(basis_->size(), basis_->size());
      RequirementType requirements;
      requirements.matrixAffordances = p_matrixAffordances;
      requirements.sols.emplace_back(displacement);
      requirements.parameter.insert({Ikarus::FEParameter::loadfactor, lambda});
      Eigen::MatrixXd matLocal;
      std::vector<GlobalIndex> dofs;
      for (auto& fe : feContainer) {
        matLocal.setZero(fe.size(),fe.size());
        fe.calculateMatrix(requirements,matLocal);
        dofs.resize(0);
        fe.globalIndices(dofs);
        for (auto i = 0; auto idi : dofs) {
          for (auto j = 0; auto idj : dofs) {
            mat(idi[0], idj[0]) += matLocal(i, j);
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

    Eigen::VectorXd& getVectorImpl(const Ikarus::VectorAffordances& p_vectorAffordances,
                                   const Eigen::VectorXd& displacement, const double& lambda) {
      vec.setZero(basis_->size());
      auto localView = basis_->localView();
      RequirementType requirements;
      requirements.vectorAffordances = p_vectorAffordances;
      requirements.sols.emplace_back(displacement);
      requirements.parameter.insert({Ikarus::FEParameter::loadfactor, lambda});
      Eigen::VectorXd vecLocal;
      std::vector<GlobalIndex> dofs;
      for (auto& fe : feContainer) {
        vecLocal.setZero(fe.size());
        dofs.resize(0);
        fe.calculateVector(requirements,vecLocal);
        fe.globalIndices(dofs);
        for (int i = 0; auto id : dofs) {
          vec(id[0]) += vecLocal(i);
          ++i;
        }
      }
      for (auto i = 0U; i < basis_->size(); ++i) {
        if (dirichletFlags->at(i)) vec[i] = 0;
      }

      return vec;
    }

    double getScalarImpl(const Ikarus::ScalarAffordances& p_scalarAffordances, const Eigen::VectorXd& displacement,
                         const double& lambda) {
      double scalar = 0.0;
      RequirementType requirements;
      requirements.scalarAffordances = p_scalarAffordances;
      requirements.sols.emplace_back(displacement);
      requirements.parameter.insert({Ikarus::FEParameter::loadfactor, lambda});

      for (auto& fe : feContainer)
        scalar += fe.calculateScalar(requirements);

      return scalar;
    }
    Basis const* basis_;
    FEContainer const& feContainer;
    std::vector<bool> const* dirichletFlags;
    Eigen::MatrixXd mat{};
    Eigen::VectorXd vec{};
  };

  template <typename Basis, typename FEContainer>  // requires Ikarus::Concepts::FlatIndexBasis<BasisEmbedded>
  class DenseFlatAssembler {
  public:
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    using GlobalIndex = typename FEContainer::value_type::GlobalIndex;
    explicit DenseFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {
      constraintsBelow_.reserve(basis_->size());
      size_t counter = 0;
      for (auto iv : std::ranges::iota_view{size_t(0), basis_->size()}) {
        constraintsBelow_.emplace_back(counter);
        if (dirichFlags[iv]) ++counter;
      }
      fixedDofs = std::ranges::count(dirichFlags, true);
    }

    Eigen::MatrixXd& getMatrix(const RequirementType& fErequirements) { return getMatrixImpl(fErequirements); }

    Eigen::VectorXd& getVector(const RequirementType& fErequirements) { return getVectorImpl(fErequirements); }

    double getScalar(const RequirementType& fErequirements) { return getScalarImpl(fErequirements); }

    Eigen::MatrixXd& getReducedMatrix(const RequirementType& fErequirements) {
      return getReducedMatrixImpl(fErequirements);
    }

    Eigen::VectorXd& getReducedVector(const RequirementType& fErequirements) {
      return getReducedVectorImpl(fErequirements);
    }

    auto createFullVector(const Eigen::VectorXd& reducedVector) {
      assert(reducedVector.size() == static_cast<Eigen::Index>(this->reducedSize())
             && "The reduced vector you passed has the wrong dimensions.");
      Eigen::Index reducedCounter = 0;
      Eigen::VectorXd fullVec(basis_->size());
      for (Eigen::Index i = 0; i < fullVec.size(); ++i) {
        if (dirichletFlags->at(i)) {
          ++reducedCounter;
          fullVec[i] = 0.0;
          continue;
        } else
          fullVec[i] = reducedVector[i - reducedCounter];
      }
      return fullVec;
    }

  private:
    Eigen::VectorXd& getVectorImpl(const RequirementType& fErequirements) {
      vec.setZero(basis_->size());
      for (auto& fe : feContainer)
        vec(fe.globalIndices()) += fe.calculateVector(fErequirements);

      return vec;
    }

    Eigen::VectorXd& getReducedVectorImpl(const RequirementType& fErequirements) {
      vecRed.setZero(reducedSize());
      int reducedCounter = 0;
      Eigen::VectorXd vecLocal;
      std::vector<GlobalIndex> dofs;
      for (auto& fe : feContainer) {
        vecLocal.setZero(fe.size());
        dofs.resize(0);
        fe.calculateVector(fErequirements,vecLocal);
        fe.globalIndices(dofs);
        assert(static_cast<long int>(dofs.size()) == vecLocal.size() && "The returned vector has wrong rowSize!");
        for (int i = 0; auto&& dofIndex : dofs) {
          if (dirichletFlags->at(dofIndex[0])) {
            ++reducedCounter;
            ++i;
            continue;
          } else
            vecRed(dofIndex[0] - constraintsBelow_[dofIndex[0]]) += vecLocal[i++];
        }
      }
      return vecRed;
    }
    Eigen::MatrixXd& getReducedMatrixImpl(const RequirementType& fErequirements) {
      matRed.setZero(reducedSize(), reducedSize());
      Eigen::MatrixXd matLocal;
      std::vector<GlobalIndex> dofs;
      for (auto& fe : feContainer) {
        matLocal.setZero(fe.size(),fe.size());
        dofs.resize(0);
        fe.calculateMatrix(fErequirements,matLocal);
        fe.globalIndices(dofs);
        assert(dofs.size() == static_cast<unsigned>(matLocal.rows()) && "The returned matrix has wrong rowSize!");
        assert(dofs.size() == static_cast<unsigned>(matLocal.cols()) && "The returned matrix has wrong colSize!");
        for (auto r = 0U; r < dofs.size(); ++r) {
          if (dirichletFlags->at(dofs[r][0])) {
            continue;
          } else {
            for (auto c = 0U; c < dofs.size(); ++c) {
              if (dirichletFlags->at(dofs[c][0])) {
                continue;
              }
              matRed(dofs[r][0] - constraintsBelow_[dofs[r][0]], dofs[c][0] - constraintsBelow_[dofs[c][0]]) += matLocal(r, c);
            }
          }
        }
      }
      return matRed;
    }

    size_t reducedSize(int i = 0) { return basis_->size() - fixedDofs; }

    size_t size(int i = 0) { return basis_->size(); }

    Eigen::MatrixXd& getMatrixImpl(const RequirementType& fErequirements) {
      mat.setZero(basis_->size(), basis_->size());
      Eigen::MatrixXd matLocal;
      std::vector<GlobalIndex> dofs;
      for (auto& fe : feContainer) {
        matLocal.setZero(fe.size(),fe.size());
        dofs.resize(0);
        fe.calculateMatrix(fErequirements,matLocal);
        fe.globalIndices(dofs);
        for (auto i = 0; auto idi : dofs) {
          for (auto j = 0; auto idj : dofs) {
            mat(idi[0], idj[0]) += matLocal(i, j);
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

    double getScalarImpl(const RequirementType& fErequirements) {
      double scalar = 0.0;

      for (auto& fe : feContainer)
        scalar += fe.calculateScalar(fErequirements);

      return scalar;
    }
    Basis const* basis_;
    FEContainer const& feContainer;
    std::vector<bool> const* dirichletFlags;
    std::vector<size_t> constraintsBelow_;
    Eigen::MatrixXd mat{};
    Eigen::MatrixXd matRed{};
    Eigen::VectorXd vec{};
    Eigen::VectorXd vecRed{};
    size_t fixedDofs{};
  };

}  // namespace Ikarus