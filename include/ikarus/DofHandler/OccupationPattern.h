//
// Created by Alex on 19.05.2021.
//

#pragma once

#include <Eigen/Sparse>
#include <memory>

template <typename GridType, typename MatrixType> class SimpleOccupationPattern {
  //    void generate()
  //    void get()

private:
  std::vector<Eigen::Triplet<double> > triplet;
  std::shared_ptr<GridType> grid;
};
