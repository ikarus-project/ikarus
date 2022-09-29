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

#include <Eigen/Core>

template <typename Derived, typename OtherDerived>
requires(std::convertible_to<Derived, Eigen::EigenBase<Derived> const&>
and std::convertible_to<OtherDerived, Eigen::EigenBase<OtherDerived> const&>)
bool
isApproxSame(Derived const& val, OtherDerived const& other, double prec) {
    if constexpr (requires {
            val.isApprox(other, prec);
            (val - other).isMuchSmallerThan(1, prec);
    })
        return val.isApprox(other, prec) or (val - other).isZero(prec);
    else if constexpr (requires { val.isApprox(other, prec); })
        return val.isApprox(other, prec);
    else  // Eigen::DiagonalMatrix branch
        return val.diagonal().isApprox(other.diagonal(), prec) or (val.diagonal() - other.diagonal()).isZero(prec);
}
