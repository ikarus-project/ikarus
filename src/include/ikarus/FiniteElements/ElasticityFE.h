//
// Created by Alex on 10.05.2021.
//

#pragma once
#include <concepts>

#include <dune/common/classname.hh>

#include <spdlog/spdlog.h>

#include <ikarus/Geometries/GeometryInterface.h>
#include <ikarus/Grids/GridEntities/GridEntitiesInterface.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>

namespace Ikarus::FiniteElements {
  template <Ikarus::Concepts::GridEntity GridEntityType, std::floating_point ct = double> class ElasticityFE {
  public:
    ElasticityFE(GridEntityType& gE) : elementGridEntity{&gE} {}

    /** \brief Type used for coordinates */
    using ctype = ct;

    /** \brief Type of the Geometry */
    using Geometry = typename GridEntityType::Geometry;

    /** \brief Dimension of the world space */
    static constexpr int coorddimension = Geometry::coorddimension;

    /** \brief Dimension of the geometry */
    static constexpr int mydimension = Geometry::mydimension;

    /** \brief Type of the Nodes */
    using NodeType = Eigen::Matrix<ctype, coorddimension, 1>;

    /** \brief Type of the ParameterSpace coordinate */
    using ParameterSpaceType = Eigen::Matrix<ctype, mydimension, 1>;

    /** \brief Type of the DofVector */
    using DofVectorType = DynArrayXi;

    /** \brief Type of the Dofs / SolutionType
     * using NodalSolutionType = Displacement<ctype,coorddimension>;*/

    /** \brief Type of the internal forces */
    using VectorType = DynVectord;

    /** \brief Type of the stiffness matrix */
    using MatrixType = DynMatrixd;

    [[nodiscard]] constexpr int dofSize() const { return GridEntityType::dimension; }

    void generateDofs() const {
//          for(auto vert: vertices(elementGridEntity))

//           dofVector.push_back(vert.addDof<Ikarus::DISPLACEMENTD_3D>())
      //        std::vector<std::shared_ptr<GenericVariableOwner>> vec;
      //        for (auto&& node : nodes)
      //            node->addVariable<Ikarus::DISPLACEMENTD_3D>();
      //
      //        for (auto&& edge : edges)
      //            edge->addVariable<Ikarus::DISPLACEMENTD_3D>();
      //
      //        for (auto&& face : faces)
      //            edge->addVariable<Ikarus::DISPLACEMENTD_3D>();

      //        this->addVariable<Ikarus::EAS>();
    }

    void initialize() { std::cout << "initialize ElasticityFE" << std::endl; }

    [[nodiscard]] std::pair<VectorType, MatrixType> calculateLocalSystem() const {
      return calculateStiffnessMatrixAndInternalForcesImpl();
    }

    [[nodiscard]] MatrixType calculateLHS() const {
      return calculateStiffnessMatrixAndInternalForcesImpl<false, true>();
    }

    [[nodiscard]] ctype getEnergy() const { return 0.0; }

    [[nodiscard]] VectorType calculateRHS() const {
      return calculateStiffnessMatrixAndInternalForcesImpl<true, false>();
    }

    template <bool internalForcesFlag = true, bool stiffnessMatrixFlag = true>
    auto calculateStiffnessMatrixAndInternalForcesImpl() const {
      if constexpr (internalForcesFlag && stiffnessMatrixFlag) {
        const VectorType Fint = VectorType::Ones(5);
        const MatrixType K = MatrixType::Ones(5, 5);
        return std::make_pair(Fint, K);
      }

      else if constexpr (internalForcesFlag && !stiffnessMatrixFlag)
        return VectorType::Ones(5);
      else if constexpr (!internalForcesFlag && stiffnessMatrixFlag)
        return MatrixType::Ones(5, 5);
      else
        static_assert(internalForcesFlag == false && stiffnessMatrixFlag == false,
                      "You asked the element: \"Don't return anything\"");
    }

    [[nodiscard]] DofVectorType getDofVector() const {
      // return localDofHandler->getDofVector();
      DofVectorType dof;
      dof.setLinSpaced(0, 5);
      //        for (auto&& node: nodes) {
      //            dof.emplace_back(node->getDofs(POSITION).ID());
      //            dof.emplace_back(node->getDofs(DIRECTOR).ID());
      //            dof.emplace_back(node->getDofs(TEMPERATURE).ID());
      //        }
      //        dof[13] = EASPARAMS;
      return dof;
    }

  private:
    GridEntityType* elementGridEntity;
  };

}  // namespace Ikarus::PhysicalElements
