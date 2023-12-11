// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
namespace Ikarus {

  namespace Impl {
    struct StressIndexPair {
      Eigen::Index row;
      Eigen::Index col;
    };

    template <size_t size>
    consteval auto createfreeVoigtIndices(const std::array<StressIndexPair, size> &fixed) {
      std::array<size_t, 6 - size> res{};
      std::array<size_t, size> voigtFixedIndices;
      std::ranges::transform(fixed, voigtFixedIndices.begin(), [](auto pair) { return toVoigt(pair.row, pair.col); });
      std::ranges::sort(voigtFixedIndices);
      std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(6)), voigtFixedIndices, res.begin());
      std::ranges::sort(res);
      return res;
    }

    template <size_t size>
    consteval auto createFixedVoigtIndices(const std::array<StressIndexPair, size> &fixed) {
      std::array<size_t, size> fixedIndices;
      std::ranges::transform(fixed, fixedIndices.begin(), [](auto pair) { return toVoigt(pair.row, pair.col); });
      std::ranges::sort(fixedIndices);
      return fixedIndices;
    }

    template <size_t size>
    constexpr size_t countDiagonalIndices(const std::array<StressIndexPair, size> &fixed) {
      size_t count = 0;
      for (auto v : fixed) {
        if (v.col == v.row) ++count;
      }
      return count;
    }
  }  // namespace Impl

  template <auto stressIndexPair, typename MaterialImpl>
  struct VanishingStress : public Material<VanishingStress<stressIndexPair, MaterialImpl>> {
    explicit VanishingStress(MaterialImpl mat, typename MaterialImpl::ScalarType p_tol = 1e-12)
        : matImpl{mat}, tol{p_tol} {}

    using Underlying = MaterialImpl;

    static constexpr auto fixedPairs                    = stressIndexPair;
    static constexpr auto freeVoigtIndices              = createfreeVoigtIndices(fixedPairs);
    static constexpr auto fixedVoigtIndices             = createFixedVoigtIndices(fixedPairs);
    static constexpr auto fixedDiagonalVoigtIndicesSize = countDiagonalIndices(fixedPairs);
    static constexpr auto freeStrains                   = freeVoigtIndices.size();
    using ScalarType                                    = typename MaterialImpl::ScalarType;
    // https://godbolt.org/z/hcs7j5rq7

    [[nodiscard]] constexpr std::string nameImpl() const noexcept {
      auto matName = matImpl.name() + "_Vanishing(";
      for (auto p : fixedPairs)
        matName += "(" + std::to_string(p.row) + std::to_string(p.col) + ")";
      matName += ")";
      return matName;
    }

    static constexpr auto strainTag          = MaterialImpl::strainTag;
    static constexpr auto stressTag          = MaterialImpl::stressTag;
    static constexpr auto tangentModuliTag   = MaterialImpl::tangentModuliTag;
    static constexpr bool energyAcceptsVoigt = MaterialImpl::energyAcceptsVoigt;
    static constexpr bool stressToVoigt      = true;
    static constexpr bool stressAcceptsVoigt = true;
    static constexpr bool moduliToVoigt      = true;
    static constexpr bool moduliAcceptsVoigt = true;

    template <typename Derived>
    ScalarType storedEnergyImpl(const Eigen::MatrixBase<Derived> &E) const {
      const auto [nonOp, Esol] = reduceStress(E);
      return matImpl.storedEnergyImpl(Esol);
    }

    template <bool voigt, typename Derived>
    auto stressesImpl(const Eigen::MatrixBase<Derived> &E) const {
      const auto [nonOp, Esol] = reduceStress(E);
      auto stressesRed         = matImpl.template stresses<MaterialImpl::strainTag, true>(Esol);

      if constexpr (voigt) {
        return removeCol(stressesRed, fixedVoigtIndices);
      } else
        return fromVoigt(stressesRed, false);
    }

    template <bool voigt, typename Derived>
    auto tangentModuliImpl(const Eigen::MatrixBase<Derived> &E) const {
      const auto [nonOp, Esol] = reduceStress(E);
      auto C                   = matImpl.template tangentModuli<MaterialImpl::strainTag, true>(Esol);
      if constexpr (voigt)
        return staticCondensation(C, fixedVoigtIndices);
      else
        return fromVoigt(C);
    }

    template <typename ScalarTypeOther>
    auto rebind() const {
      auto reboundMatImpl = matImpl.template rebind<ScalarTypeOther>();
      return VanishingStress<stressIndexPair, decltype(reboundMatImpl)>(reboundMatImpl, tol);
    }

  private:
    template <typename Derived>
    decltype(auto) maybeFromVoigt(const Eigen::MatrixBase<Derived> &E) const {
      if constexpr (Concepts::EigenVector<Derived>) {  // receiving vector means Voigt notation
        return fromVoigt(E.derived(), true);
      } else
        return E.derived();
    }

    template <typename Derived>
    void initUnknownStrains(Eigen::MatrixBase<Derived> &E) const {
      for (size_t i = 0; i < fixedPairs.size(); ++i) {
        ScalarType initialVal = E(fixedPairs[i].row, fixedPairs[i].col);
        if constexpr (strainTag == StrainTags::deformationGradient or strainTag == StrainTags::rightCauchyGreenTensor) {
          if (Dune::FloatCmp::eq(initialVal, ScalarType(0.0)) and (fixedPairs[i].row == fixedPairs[i].col))
            initialVal = ScalarType(1.0);
        }
        if (fixedPairs[i].row != fixedPairs[i].col) initialVal = ScalarType(0.0);
        E(fixedPairs[i].row, fixedPairs[i].col) = E(fixedPairs[i].col, fixedPairs[i].row) = initialVal;
      }
    }

    template <typename Derived>
    auto reduceStress(const Eigen::MatrixBase<Derived> &p_Eraw) const {
      auto E = maybeFromVoigt(p_Eraw);
      initUnknownStrains(E);

      std::array<size_t, fixedDiagonalVoigtIndicesSize> fixedDiagonalVoigtIndices;
      for (size_t ri = 0; auto i : fixedVoigtIndices) {
        auto indexPair = fromVoigt(i);
        if (indexPair[0] == indexPair[1]) fixedDiagonalVoigtIndices[ri++] = i;
      }
      auto f = [&](auto &) {
        auto S = matImpl.template stresses<MaterialImpl::strainTag, true>(E);
        return S(fixedDiagonalVoigtIndices).eval();
      };
      auto df = [&](auto &) {
        auto moduli = (matImpl.template tangentModuli<MaterialImpl::strainTag, true>(E)).eval();
        return (moduli(fixedDiagonalVoigtIndices, fixedDiagonalVoigtIndices) / MaterialImpl::derivativeFactor).eval();
      };

      auto Er    = E(fixedDiagonalVoigtIndices, fixedDiagonalVoigtIndices).eval().template cast<ScalarType>();
      auto nonOp = Ikarus::NonLinearOperator(functions(f, df), parameter(Er));
      auto nr    = Ikarus::makeNewtonRaphson(
             nonOp, [&](auto &r, auto &A) { return (A.inverse() * r).eval(); },
             [&](auto    &/* Ex33 */, auto &Ecomps) {
            for (int ri = 0; auto i : fixedDiagonalVoigtIndices) {
              auto indexPair = fromVoigt(i);
              E(indexPair[0], indexPair[1]) += Ecomps(ri++);
            }
             });
      nr->setup({.tol = tol, .maxIter = 100});
      if (!static_cast<bool>(nr->solve()))
        DUNE_THROW(Dune::MathError, "The stress reduction of material " << nameImpl() << " was unsuccessful\n"
                                                                        << "The strains are\n"
                                                                        << E << "\n The stresses are\n"
                                                                        << f(Er));
      return std::make_pair(nonOp, E);
    }

    MaterialImpl matImpl;
    double tol{};
  };

  template <Impl::StressIndexPair... stressIndexPair, typename MaterialImpl>
  auto makeVanishingStress(MaterialImpl mat, typename MaterialImpl::ScalarType p_tol = 1e-12) {
    return VanishingStress<std::to_array({stressIndexPair...}), MaterialImpl>(mat, p_tol);
  }

  template <typename MaterialImpl>
  auto planeStress(const MaterialImpl &mat, typename MaterialImpl::ScalarType p_tol = 1e-8) {
    return makeVanishingStress<Impl::StressIndexPair{2, 1}, Impl::StressIndexPair{2, 0}, Impl::StressIndexPair{2, 2}>(
        mat, p_tol);
  }

  template <typename MaterialImpl>
  auto shellMaterial(const MaterialImpl &mat, typename MaterialImpl::ScalarType p_tol = 1e-8) {
    return makeVanishingStress<Impl::StressIndexPair{2, 2}>(mat, p_tol);
  }

  template <typename MaterialImpl>
  auto beamMaterial(const MaterialImpl &mat, typename MaterialImpl::ScalarType p_tol = 1e-8) {
    return makeVanishingStress<Impl::StressIndexPair{1, 1}, Impl::StressIndexPair{2, 2}>(mat, p_tol);
  }
}  // namespace Ikarus
