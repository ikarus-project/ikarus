//
// Created by alex on 12/20/21.
//

#pragma once
#include <functional>
#include <map>
#include <ranges>

#include <Eigen/Core>

namespace Ikarus {
  template <typename FEManager>
  class DirichletConditionManager {
  public:
    static constexpr int gridDim = 2;
    //    using SpaceFunctionType      = std::function<double(Eigen::Vector<double, gridDim>&)>;
    using TimeFunctionType = std::function<double(double&)>;
    explicit DirichletConditionManager(const FEManager& feManager) : feManager_{&feManager} {}

    template <class GridEntity>
    void addConstraint(
        const GridEntity& gridEntity, const int localVariableIndex,
        //        const SpaceFunctionType spaceFunctionType = []([[maybe_unused]] auto&& v) { return 1.0; },
        const TimeFunctionType timeFunctionType = []([[maybe_unused]] auto&& v) { return 1.0; }) {
//
//      auto gridEntityIndices = feManager_->dofIndicesOfEntity(gridEntity);
//      std::cout<<"gridEntityIndices.transpose()"<<std::endl;
//      std::cout<<gridEntityIndices.transpose()<<std::endl;
//      std::cout<<"localVariableIndex"<<std::endl;
//      std::cout<<localVariableIndex<<std::endl;
//      std::cout<<"gridEntityIndices[localVariableIndex]"<<std::endl;
//      std::cout<<gridEntityIndices[localVariableIndex]<<std::endl;
//      assert(
//          localVariableIndex < gridEntityIndices.size()
//          && "The index of the degree of freedom you want to constrain is larger then the underlying variable dof "
//             "count.");
//
//
//      hasDirichletBoundaryCondition_.insert(
//          {gridEntityIndices[localVariableIndex], {gridEntityIndices[localVariableIndex], timeFunctionType}});
//      std::cout<<"========="<<std::endl;
      isFinalized_ = false;
    }

    [[nodiscard]] bool isConstrained(size_t i) const { return hasDirichletBoundaryCondition_.contains(i); }

    [[nodiscard]] size_t constraintsBelow(size_t i) const {
      assert(isFinalized_ && "You have to call the finalize method before you can use constraintsBelow");
      return constraintsBelow_[i];
    }

    void finalize() {
      if (not isFinalized_) {
        constraintsBelow_.reserve(feManager_->numberOfDegreesOfFreedom());
        for (auto iv : std::ranges::iota_view{size_t(0), feManager_->numberOfDegreesOfFreedom()}) {
          constraintsBelow_.emplace_back(
              std::distance(hasDirichletBoundaryCondition_.begin(), hasDirichletBoundaryCondition_.lower_bound(iv)));
        }
      }
      isFinalized_ = true;
    }

    auto getValue(size_t i) const {
      assert(isConstrained(i));
      return hasDirichletBoundaryCondition_[i];
    }

    auto numberOfReducedDegreesOfFreedom() const {
      return feManager_->numberOfDegreesOfFreedom() - hasDirichletBoundaryCondition_.size();
    }

    auto constrainedIndices() { return std::views::keys(hasDirichletBoundaryCondition_); }

    auto freeIndices() {
      return std::ranges::iota_view{size_t(0), feManager_->numberOfDegreesOfFreedom()}
             | std::views::filter([&](auto&& i) { return !isConstrained(i); });
    }

    auto viewAsFullVector(const Eigen::VectorXd& v) {
      return std::ranges::iota_view{std::size_t{0}, feManager_->numberOfDegreesOfFreedom()}
             | std::views::transform([this, &v](int i) mutable {
                 if (this->isConstrained(i)) {
                   return 0.0;
                 } else
                   return v[i - this->constraintsBelow(i)];
               });
    }

  private:
    std::vector<size_t> constraintsBelow_;
    bool isFinalized_{false};
    struct DirichletConstrainedDof {
      size_t index{};
      //      SpaceFunctionType spaceFunctionType;
      TimeFunctionType timeFunctionType;
    };
    std::map<size_t, DirichletConstrainedDof> hasDirichletBoundaryCondition_;
    FEManager const* feManager_;
  };
}  // namespace Ikarus