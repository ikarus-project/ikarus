
#pragma once
#include <concepts>

#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
namespace Ikarus::FiniteElements {

  template <typename GridElementEntityType, typename IndexSetType, std::floating_point ct = double>
  class ForceLoad : public FEVertexDisplacement<GridElementEntityType, IndexSetType> {
  public:
    using Traits = FETraits<GridElementEntityType>;
    using Base   = FEVertexDisplacement<GridElementEntityType, IndexSetType>;
    using Base::getEntityVariableTuple;
    ForceLoad(GridElementEntityType &gE, const IndexSetType &indexSet)
        : Base(gE, indexSet), elementGridEntity{&gE}, indexSet_{&indexSet} {}

    [[nodiscard]] typename Traits::VectorType calculateVector(const typename Traits::FERequirementType &req) const {
      return calculateVectorImpl(req);
    }
    [[nodiscard]] typename Traits::MatrixType calculateMatrix(
        [[maybe_unused]] const typename Traits::FERequirementType &req) const {
      return typename Traits::MatrixType(this->dofSize(), this->dofSize());
    }

    [[nodiscard]] typename Traits::VectorType calculateVectorImpl(
        [[maybe_unused]] const typename Traits::FERequirementType &req) const {
      assert(req.parameter.contains(FEParameter::loadfactor));
      const auto rule = Dune::QuadratureRules<double, Traits::mydim>::rule(duneType(elementGridEntity->type()), 2);
      typename Traits::VectorType Fext(this->dofSize());
      Fext.setZero();

      auto f = [&]([[maybe_unused]] auto &gpPos) {
        Eigen::Vector<double, Traits::dimension> feval;
        feval.setZero();
        auto globalCoords = elementGridEntity->geometry().global(gpPos);
        feval[1]          = 1.0 * getValue(req.parameter.at(FEParameter::loadfactor))[0] * 1;  // globalCoords [0]

        return feval;
      };
      for (auto &gp : rule) {
        const auto N = Ikarus::LagrangeCube<double, Traits::mydim, 1>::evaluateFunction(toEigenVector(gp.position()));

        for (int vertexCounter = 0; auto Ni : N) {
          Fext.template segment<Traits::dimension>(vertexCounter) -= Ni * f(gp.position());
          vertexCounter += Traits::dimension;
        }
      }
      return Fext;
    }

  private:
    GridElementEntityType const *const elementGridEntity;
    IndexSetType const *const indexSet_;
  };

}  // namespace Ikarus::FiniteElements