
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
    using GlobalCoordinates = typename Traits::GlobalCoordinates;
    template<typename SpaceFunction,typename TimeFunctionReturnType=double,typename TimeFunctionParameterType=double>
    ForceLoad(GridElementEntityType &gE, const IndexSetType &indexSet,SpaceFunction&& spaceFunction,std::function<TimeFunctionReturnType(TimeFunctionParameterType)> timeFunction = []( double t)-> double{return t;})
        : Base(gE, indexSet), elementGridEntity{&gE}, indexSet_{&indexSet},spaceFunction_{spaceFunction},timeFunction_{timeFunction} {}

    [[nodiscard]] typename Traits::VectorType calculateVector(const typename Traits::FERequirementType &req) const {
      return calculateVectorImpl(req);
    }
    [[nodiscard]] typename Traits::MatrixType calculateMatrix(
        [[maybe_unused]] const typename Traits::FERequirementType &req) const {
      return Traits::MatrixType::Zero(this->dofSize(), this->dofSize());
    }

    [[nodiscard]] typename Traits::VectorType calculateVectorImpl(
        [[maybe_unused]] const typename Traits::FERequirementType &req) const {
      assert(req.parameter.contains(FEParameter::time));
      const auto rule = Dune::QuadratureRules<double, Traits::mydim>::rule(duneType(elementGridEntity->type()), 2);
      typename Traits::VectorType Fext(this->dofSize());
      Fext.setZero();
      const auto ft= timeFunction_(getValue(req.parameter.at(FEParameter::time))[0]);
      for (auto &gp : rule) {
        const auto N = Ikarus::LagrangeCube<double, Traits::mydim, 1>::evaluateFunction(toEigenVector(gp.position()));

        for (int vertexCounter = 0; auto Ni : N) {
          Fext.template segment<Traits::dimension>(vertexCounter) -= Ni * spaceFunction_(toEigenVector(elementGridEntity->geometry().global(gp.position())))*ft;
          vertexCounter += Traits::dimension;
        }
      }
      return Fext;
    }

  private:
    GridElementEntityType const *const elementGridEntity;
    IndexSetType const *const indexSet_;
    std::function<GlobalCoordinates(const GlobalCoordinates&)> spaceFunction_;
    std::function<double(double)> timeFunction_;
  };

}  // namespace Ikarus::FiniteElements