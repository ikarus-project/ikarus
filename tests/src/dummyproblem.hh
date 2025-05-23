// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <config.h>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>

template <typename G = Dune::UGGrid<2>, bool useYASP = false>
struct DummyProblem
{
  using Grid     = G;
  using GridView = typename Grid::LeafGridView;

  using PreBasis = Dune::Functions::PowerPreBasis<Dune::Functions::BasisFactory::FlatInterleaved,
                                                  Dune::Functions::LagrangePreBasis<GridView, 1>, 2ul>;
  using Basis    = Ikarus::BasisHandler<PreBasis>;

  using Material = typename Ikarus::Materials::VanishingStress<std::array<Ikarus::Materials::MatrixIndexPair, 3ul>{
                                                                   {Ikarus::Materials::MatrixIndexPair{2ul, 1ul},
                                                                    Ikarus::Materials::MatrixIndexPair{2ul, 0ul},
                                                                    Ikarus::Materials::MatrixIndexPair{2ul, 2ul}}
  },
                                                               Ikarus::Materials::LinearElasticityT<double>>;
  using LinearElastic =
      Ikarus::FE<Ikarus::PreFE<Basis>, Ikarus::LinearElasticPre<Material>::Skill, Ikarus::VolumeLoadPre<2>::Skill>;

  using SparseAssmblerT =
      Ikarus::SparseFlatAssembler<std::vector<LinearElastic>&, Ikarus::DirichletValues<typename Basis::FlatBasis>>;

  // YASPGrid needs an int, structuresgridfactory an unsigned int
  explicit DummyProblem(
      const std::array<std::conditional_t<useYASP, int, unsigned int>, 2>& elementsPerDirection = {10, 10})
      : grid_([&]() {
          constexpr double Lx                     = 4.0;
          constexpr double Ly                     = 4.0;
          const Dune::FieldVector<double, 2> bbox = {Lx, Ly};

          if constexpr (not useYASP)
            return Dune::StructuredGridFactory<Grid>::createCubeGrid({0, 0}, bbox, elementsPerDirection);
          else
            return make_unique<Grid>(bbox, elementsPerDirection);
        }()),
        gridView_([&]() { return grid_->leafGridView(); }()),
        basis_([&]() {
          using namespace Dune::Functions::BasisFactory;

          return Ikarus::makeBasis(gridView_, power<2>(lagrange<1>(), FlatInterleaved{}));
        }()),
        dirichletValues_([&]() { return Ikarus::DirichletValues(basis_.flat()); }()),
        fes_([&]() {
          dirichletValues_.fixBoundaryDOFs(
              [&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
                if (std::abs(intersection.geometry().center()[1]) < 1e-8)
                  dirichletFlags[localView.index(localIndex)] = true;
              });
          auto vL     = []([[maybe_unused]] auto& globalCoord, auto& lamb) { return Eigen::Vector2d{0, -1}; };
          auto linMat = Ikarus::Materials::LinearElasticity(
              Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 100, .nu = 0.2}));
          auto skills_ =
              Ikarus::skills(Ikarus::linearElastic(Ikarus::Materials::planeStress(linMat)), Ikarus::volumeLoad<2>(vL));
          std::vector<LinearElastic> fes;

          for (auto&& element : elements(gridView_)) {
            fes.emplace_back(Ikarus::makeFE(basis_, skills_));
            fes.back().bind(element);
          }
          return std::move(fes);
        }()),
        sparseAssembler_{std::make_shared<SparseAssmblerT>(fes_, dirichletValues_)},
        requirement_(typename LinearElastic::Requirement(basis_))

  {
    requirement_.parameter() = 1.0;
    sparseAssembler_->bind(requirement_);
    auto f = Ikarus::DifferentiableFunctionFactory::op(
        sparseAssembler_,
        Ikarus::AffordanceCollection(Ikarus::VectorAffordance::forces, Ikarus::MatrixAffordance::stiffness));

    const auto& K    = derivative(f)(requirement_);
    const auto& Fext = f(requirement_);

    auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
    linSolver.compute(K);
    linSolver.solve(requirement_.globalSolution(), -Fext);
    sparseAssembler_->bind(requirement_); // the requirement has changed therefore, we have to bind again
    sparseAssembler_->bind(Ikarus::DBCOption::Full);
  }

  const auto& grid() { return *grid_; }
  auto& gridView() { return gridView_; }
  auto& basis() { return basis_; }
  auto& requirement() { return requirement_; }
  auto& sparseAssembler() { return sparseAssembler_; }
  auto& dirichletValues() { return dirichletValues_; }
  auto& finiteElements() { return fes_; }

private:
  std::unique_ptr<Grid> grid_;
  GridView gridView_;
  Basis basis_;
  Ikarus::DirichletValues<typename Basis::FlatBasis> dirichletValues_;
  std::vector<LinearElastic> fes_;
  std::shared_ptr<SparseAssmblerT> sparseAssembler_;
  typename LinearElastic::Requirement requirement_{};
};
