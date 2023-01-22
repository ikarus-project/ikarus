// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <dune/python/pybind11/pybind11.h>

namespace Ikarus::Python{

  template< class K, int m, int n >
  inline static void registertoVoigt ( pybind11::handle scope )
  {
    typedef Dune::FieldMatrix< K, m, n > FM;

    std::string fname = className<K>();
    auto entry = insertClass<FM>( scope, "FieldMatrix_"+fname+"_"+std::to_string(m)+"_"+std::to_string(n), pybind11::buffer_protocol(),
                                 GenerateTypeName("Dune::FieldMatrix",Dune::MetaType<K>(),m,n), IncludeFiles{"dune/common/fmatrix.hh"}
    );
    if (!entry.second)
      return;
    registerFieldMatrix( scope, entry.first );
  }
} // namespace Python

}
