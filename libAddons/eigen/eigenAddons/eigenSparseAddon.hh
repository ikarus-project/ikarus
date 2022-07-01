

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

inline Index getLinearIndex(Index row, Index col) {
  eigen_assert(row >= 0 && row < rows() && col >= 0 && col < cols());

  const Index outer = IsRowMajor ? row : col;
  const Index inner = IsRowMajor ? col : row;

  Index start = m_outerIndex[outer];
  Index end   = m_innerNonZeros ? m_outerIndex[outer] + m_innerNonZeros[outer] : m_outerIndex[outer + 1];
  eigen_assert(end >= start && "you probably called coeffRef on a non finalized matrix");

  Index p = m_data.searchLowerIndex(start, end - 1, StorageIndex(inner));
  return p;
}
