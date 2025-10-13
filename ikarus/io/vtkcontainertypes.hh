// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file vtkcontainertypes.hh
 * \brief Implementation of container types for the vtk writer.
 * \ingroup io
 *
 */

#pragma once
#include <array>
#include <tuple>

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

namespace Ikarus::Vtk::Impl {

template <typename Container>
constexpr auto sizeOfContainer = []() {
  if constexpr (requires { std::tuple_size<Container>::value; })
    return std::tuple_size<Container>::value;
  else
    return 1ul;
}();

template <class PreBasis>
struct ResultContainerPre
{
  using type = double;
};

template <class Basis>
struct ResultContainer
{
  using type = ResultContainerPre<typename Basis::PreBasis>::type;
};

// specialization for composite basis with where the subspace basis are all either power basis or scalar basis
template <class IMS, class... SPB>
struct ResultContainerPre<Dune::Functions::CompositePreBasis<IMS, SPB...>>
{
  using SubPreBases = std::tuple<SPB...>;
  using FirstBasis  = std::tuple_element_t<0, SubPreBases>;

  static constexpr size_t size = sizeof...(SPB);
  using type                   = std::array<typename ResultContainerPre<FirstBasis>::type, size>;
};

// Subspace basis can either be from a composite or a power basis
template <class PreBasis, class PP>
struct ResultContainerSSB
{
  using type = double;
};

template <class PP, class IMS, class... SPB>
struct ResultContainerSSB<Dune::Functions::CompositePreBasis<IMS, SPB...>, PP>
{
  using SubPreBases = std::tuple<SPB...>;
  using Basis       = std::tuple_element_t<PP{}.template get<0>(), SubPreBases>;

  using type = ResultContainerPre<Basis>::type;
};

// Specialization for subspace basis
template <class RB, class PP>
struct ResultContainer<Dune::Functions::SubspaceBasis<RB, PP>>
{
  using RBPre = RB::PreBasis;

  using type = ResultContainerSSB<RBPre, PP>::type;
};

template <class Basis>
using ResultContainer_t = typename ResultContainer<Basis>::type;

// specialization for power basis
template <class IMS, class SPB, std::size_t size>
struct ResultContainerPre<Dune::Functions::PowerPreBasis<IMS, SPB, size>>
{
  using type = std::array<typename ResultContainerPre<SPB>::type, size>;
};

} // namespace Ikarus::Vtk::Impl