//
// Created by Alex on 26.05.2021.
//

#pragma once

template <Grid GridType class DefaultDofHandler {
  DefaultDofHandler(shared_ptr<GridType> gridInput) : grid{gridInput} {}

private:
  std::shared_ptr<GridType> grid;
};
