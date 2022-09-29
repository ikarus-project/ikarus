/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#pragma once
#include <ikarus/localFunctions/meta.hh>
namespace Ikarus {

  template <typename LFImpl>
  class ClonableLocalFunction {
  public:
    LFImpl clone() const { return LFImpl(underlying().basis(), underlying().coeffs, typename LFImpl::Ids()); }

    template <typename OtherType, size_t ID = 0>
    auto rebindClone(OtherType&& , [[maybe_unused]] Dune::index_constant<ID>&& id = Dune::index_constant<0>()) const {
      if constexpr (LFImpl::Ids::value == ID)
        return typename LFImpl::template Rebind<OtherType>::other(
            underlying().basis(), convertUnderlying<OtherType>(underlying().coeffs), typename LFImpl::Ids());
      else
        return clone();
    }

  private:
    constexpr LFImpl const& underlying() const  // CRTP
    {
      return static_cast<LFImpl const&>(*this);
    }
  };

}  // namespace Ikarus
