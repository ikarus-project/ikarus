//
// Created by Alex on 19.05.2021.
//

#pragma once

namespace Ikarus {
  template <typename ResidualAssemblerType, typename GridType, typename DofHandlerType>
  concept ResidualAssembler = requires(ResidualAssemblerType assembler&& GridType grid) {
    typename ResidualAssemblerType::ResidualType;
    { fe.calculateLHS(grid) } -> std::same_as<typename FEConceptType::MatrixType>;
  };

};  // namespace Ikarus
