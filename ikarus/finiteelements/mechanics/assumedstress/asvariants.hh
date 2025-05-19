// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file asvariants.hh
 * \brief Definition of the AssumedStress variants.
 * \ingroup mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/assumedstress/asfunctions.hh>
#include <ikarus/finiteelements/mechanics/assumedstress/asvariants/linearandpk2stress.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::PS {
namespace Impl {
  template <typename GEO, int myDim, typename ASFunction>
  struct AssumedStressType;
  template <typename GEO>
  struct AssumedStressType<GEO, 3, PS::LinearStress>
  {
    using type = std::variant<S18<GEO>, S24<GEO>, S30<GEO>>;
  };
  template <typename GEO>
  struct AssumedStressType<GEO, 2, PS::LinearStress>
  {
    using type = std::variant<S5<GEO>>;
  };
  template <typename GEO, int myDim>
  struct AssumedStressType<GEO, myDim, PS::PK2Stress>
  {
    using type = AssumedStressType<GEO, myDim, PS::LinearStress>::type;
  };
} // namespace Impl

/**
 * \brief Wrapper around the AssumedStress variant, contains helper functions
 * \tparam GEO the Geometry type of the udderlying grid element
 * \tparam ASFunction the assumed stress function
 */
template <typename ASFunction, typename GEO>
struct AssumedStressVariant
{
  static constexpr int myDim = GEO::mydimension;
  using Variant              = Impl::AssumedStressType<GEO, myDim, ASFunction>::type;
  static_assert((myDim == 2) or (myDim == 3), "AssumedStress variants are only available for 2D and 3D elements.");

  /**
   * \brief Executes a function F on the variant
   * \tparam F the Type of the function F
   * \param f the function, which takes the as element as an argument
   */
  template <typename F>
  void operator()(F&& f) const {
    std::visit([&]<typename AssumedStressT>(const AssumedStressT& asFunction) { f(asFunction); }, var_);
  }

  /**
   * \brief A helper function to get the number of AssumedStress parameters.
   * \return Number of AssumedStress parameters.
   */
  auto numberOfInternalVariables() const {
    return std::visit(
        []<typename AssumedStressT>(const AssumedStressT&) -> int { return AssumedStressT::assumedStressSize; }, var_);
  }

  void setAssumedStressType(int numberOfInternalVariables) {
    numberOfAssumedStressParameters_ = numberOfInternalVariables;
    if (geometry_)
      createAssumedStressType();
  }
  void bind(const GEO& geometry) {
    geometry_ = std::make_optional<GEO>(geometry);
    createAssumedStressType();
  }

private:
  void createAssumedStressType() {
    const std::string& errorMessage =
        "The given AssumedStress parameters are not available for the " + std::to_string(myDim) + "D case.";

    if constexpr (std::same_as<ASFunction, PS::LinearStress> or std::same_as<ASFunction, PS::PK2Stress>) {
      if constexpr (myDim == 2) {
        switch (numberOfAssumedStressParameters_) {
          case 5:
            var_ = S5(geometry_.value());
            break;
          default:
            DUNE_THROW(Dune::NotImplemented, errorMessage);
        }
      } else if constexpr (myDim == 3) {
        switch (numberOfAssumedStressParameters_) {
          case 18:
            var_ = S18(geometry_.value());
            break;
          case 24:
            var_ = S24(geometry_.value());
            break;
          case 30:
            var_ = S30(geometry_.value());
            break;
          default:
            DUNE_THROW(Dune::NotImplemented, errorMessage);
        }
      }
    }
  }
  std::optional<GEO> geometry_;
  Variant var_;
  int numberOfAssumedStressParameters_;
};
} // namespace Ikarus::PS
