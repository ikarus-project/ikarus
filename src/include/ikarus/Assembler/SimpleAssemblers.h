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

  template <typename Basis, typename FEContainer>  // requires Ikarus::Concepts::FlatIndexBasis<Basis>
  class SparseFlatAssembler {
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    //  using ScalarAssembler = ScalarAssembler<Basis, FEContainer>;
    //  using VectorAssembler = VectorFlatAssembler<Basis, FEContainer>;

  public:
    SparseFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {
      constraintsBelow_.reserve(basis_->size());
      size_t counter = 0;
      for (auto iv : std::ranges::iota_view{size_t(0), basis_->size()}) {
        if(dirichFlags[iv])
          ++counter;
        constraintsBelow_.emplace_back(counter--);
      }

      fixedDofs = std::ranges::count(dirichFlags,true);
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
        A = fe.calculateMatrix(fErequirements);
        assert(std::sqrt(elementLinearIndices[elementIndex].size()) == A.rows() && "The returned matrix has wrong rowSize!");
        assert(std::sqrt(elementLinearIndices[elementIndex].size()) == A.cols() && "The returned matrix has wrong colSize!");
        for (Eigen::Index linearIndex = 0; double matrixEntry : A.reshaped())
          spMat.coeffs()(elementLinearIndices[elementIndex][linearIndex++]) += matrixEntry;
        ++elementIndex;
      }
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) spMat.col(i)*=0;
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) spMat.row(i)*=0;
      for (auto i = 0U; i < basis_->size(); ++i)
        if (dirichletFlags->at(i)) spMat.diagonal()[i] = 1;
      return spMat;
    }

    Eigen::SparseMatrix<double>& getReducedMatrixImpl(const RequirementType& fErequirements) {
      if (!isReducedOccupationPatternCreated) createReducedOccupationPattern();
      if (!arelinearReducedDofsPerElementCreated) createlinearDofsPerElementReduced();
      spMatReduced.coeffs().setZero();
      Eigen::MatrixXd A;
      for (size_t elementIndex = 0; const auto& fe : feContainer) {
        A               = fe.calculateMatrix(fErequirements);
        const auto dofs = fe.globalIndices();
        assert(dofs.size() == A.rows() && "The returned matrix has wrong rowSize!");
        assert(dofs.size() == A.cols() && "The returned matrix has wrong colSize!");
        Eigen::Index linearIndex = 0;
        for (int r = 0; r < dofs.size(); ++r) {
          if (dirichletFlags[dofs[r]])
            continue;
          else {
            for (int c = 0; c < dofs.size(); ++c) {
              if (dirichletFlags[dofs[c]]) continue;
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
      for (auto&& fe : feContainer)
        for (auto idi : fe.globalIndices())
          for (auto idj : fe.globalIndices())
            vectorOfTriples.emplace_back(idi, idj, 0.0);

      spMat.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
      isOccupationPatternCreated = true;
    }

    // https://stackoverflow.com/questions/59192659/efficiently-use-eigen-for-repeated-sparse-matrix-assembly-in-nonlinear-finite-el
    void createReducedOccupationPattern() {
      spMatReduced.resize(reducedSize(), reducedSize());
      std::vector<Eigen::Triplet<double>> vectorOfTriples;
      const int estimateOfConnectivity = 8;
      using std::size;

      vectorOfTriples.reserve(estimateOfConnectivity * basis_->GridView()->size(GridView::dimension));
      for (auto& fe : feContainer) {
        const auto dofs = fe.globalIndices();
        for (int r = 0; r < dofs.size(); ++r) {
          if (dirichletFlags[dofs[r]])
            continue;
          else {
            for (int c = 0; c < dofs.size(); ++c) {
              if (dirichletFlags[dofs[c]]) continue;
              vectorOfTriples.emplace_back(dofs[r] - constraintsBelow_[dofs[r]],
                                           dofs[c] - constraintsBelow_[dofs[c]], 0.0);
            }
          }
        }
      }

      spMatReduced.setFromTriplets(vectorOfTriples.begin(), vectorOfTriples.end());
      isReducedOccupationPatternCreated = true;
    }

    // This function save the indices of each element in the underlying vector which stores the sparse matrix entries
    void createlinearDofsPerElement() {
      for (auto&& dofsOfElement : dofsOfElements()) {
        elementLinearIndices.emplace_back(Dune::Power<2>::eval(dofsOfElement.size()));
        for (Eigen::Index linearIndexOfElement = 0; auto&& c : dofsOfElement)
          for (auto&& r : dofsOfElement)
            elementLinearIndices.back()[linearIndexOfElement++] = spMat.getLinearIndex(r, c);
      }
      arelinearDofsPerElementCreated = true;
    }

    // This function save the indices of each element in the underlying vector which stores the sparse matrix entries
    void createlinearDofsPerElementReduced() {
      for (auto&& dofs : dofsOfElements()) {
        elementLinearReducedIndices.emplace_back();
        for (int r = 0; r < dofs.size(); ++r) {
          if (dirichletFlags[dofs[r]]) continue;
          for (int c = 0; c < dofs.size(); ++c) {
            if (dirichletFlags[dofs[c]]) continue;
            elementLinearReducedIndices.back().push_back(
                spMatReduced.getLinearIndex(dofs[r] - constraintsBelow_(dofs[r]),
                                            dofs[c] - constraintsBelow_(dofs[c])));
          }
        }
      }
      arelinearReducedDofsPerElementCreated = true;
    }

    size_t reducedSize(int i = 0) { return basis_->size() - fixedDofs; }

    size_t size(int i = 0) { return basis_->size(); }

  private:
    auto dofsOfElements()
    {
      return feContainer | std::views::transform([](auto&& fe){ return fe.globalIndices();});
    }


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

  public:
    VectorFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {
      constraintsBelow_.reserve(basis->size());
      size_t counter = 0;
      for (auto iv : std::ranges::iota_view{size_t(0), basis->size()}) {
        if(dirichFlags[iv])
          ++counter;
        constraintsBelow_.emplace_back(counter--);
      }
      fixedDofs = std::ranges::count(dirichFlags,true);
    }

    Eigen::VectorXd& getVector(const RequirementType& fErequirements) {
      return getVectorImpl(fErequirements);
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

    size_t reducedSize(int i = 0) { return basis_->size() - fixedDofs; }

    size_t size(int i = 0) { return basis_->size(); }

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
      for (auto& fe : feContainer) {
        const auto f = fe.calculateVector(fErequirements);
        const auto dofs = f.globalIndices();
        assert(dofs.size() == f.size() && "The returned vector has wrong rowSize!");
        for (int i = 0; auto&& dofIndex : dofs) {
          if (dirichletFlags[dofIndex]) {
            ++reducedCounter;
            ++i;
            continue;
          } else
            vecRed(dofIndex - constraintsBelow_[dofIndex]) += f[i++];
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

    double& getScalar(const RequirementType& fErequirements) {
      return getScalarImpl(fErequirements);
    }

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

  template <typename Basis, typename FEContainer>  // requires Ikarus::Concepts::FlatIndexBasis<Basis>
  class DenseFlatSimpleAssembler {
  public:
    using RequirementType = typename FEContainer::value_type::FERequirementType;
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
      for (auto& fe : feContainer) {
        auto matLoc        = fe.calculateMatrix(requirements);
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

    Eigen::VectorXd& getVectorImpl(const Ikarus::VectorAffordances& p_vectorAffordances,
                                   const Eigen::VectorXd& displacement, const double& lambda) {
      vec.setZero(basis_->size());
      auto localView = basis_->localView();
      RequirementType requirements;
      requirements.vectorAffordances = p_vectorAffordances;
      requirements.sols.emplace_back(displacement);
      requirements.parameter.insert({Ikarus::FEParameter::loadfactor, lambda});
      for (auto& fe : feContainer) {
        //      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)> fe(localView, 1000,
        //      0.3);

        auto vecLocal = fe.calculateVector(requirements);
        for (int i = 0; auto id : fe.globalIndices()) {
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
      vec.setZero(basis_->size());
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

  template <typename Basis, typename FEContainer>  // requires Ikarus::Concepts::FlatIndexBasis<Basis>
  class DenseFlatAssembler {
  public:
    using RequirementType = typename FEContainer::value_type::FERequirementType;
    explicit DenseFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
        : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {}

    Eigen::MatrixXd& getMatrix(const RequirementType& fErequirements) {
      return getMatrixImpl(fErequirements);
    }

    Eigen::VectorXd& getVector(const RequirementType& fErequirements) {
      return getVectorImpl(fErequirements);
    }

    double getScalar(const RequirementType& fErequirements) {
      return getScalarImpl(fErequirements);
    }

  private:
    Eigen::MatrixXd& getMatrixImpl(const RequirementType& fErequirements) {
      mat.setZero(basis_->size(), basis_->size());
      for (auto& fe : feContainer) {
        auto matLoc        = fe.calculateMatrix(fErequirements);
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

    Eigen::VectorXd& getVectorImpl(const RequirementType& fErequirements) {
      vec.setZero(basis_->size());
      for (auto& fe : feContainer) {

        auto vecLocal = fe.calculateVector(fErequirements);
        for (int i = 0; auto id : fe.globalIndices()) {
          vec(id[0]) += vecLocal(i);
          ++i;
        }
      }
      for (auto i = 0U; i < basis_->size(); ++i) {
        if (dirichletFlags->at(i)) vec[i] = 0;
      }

      return vec;
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
    Eigen::MatrixXd mat{};
    Eigen::VectorXd vec{};
  };

}  // namespace Ikarus