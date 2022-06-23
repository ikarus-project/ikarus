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
#include <ostream>

namespace Ikarus::Concepts {
  template <typename ManifoldType>
  concept Manifold = requires(ManifoldType var, typename ManifoldType::CorrectionType correction, std::ostream& s,
                              typename ManifoldType::CoordinateType value, int i) {
    typename ManifoldType::ctype;
    ManifoldType::valueSize;
    ManifoldType::correctionSize;
    typename ManifoldType::CoordinateType;
    typename ManifoldType::CorrectionType;
    { var.getValue() } -> std::convertible_to<typename ManifoldType::CoordinateType>;
    { var.setValue(value) } -> std::same_as<void>;
    { var += correction } -> std::same_as<void>;
    { s << var } -> std::same_as<std::ostream&>;
  };
}  // namespace Ikarus::Concepts