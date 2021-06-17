//
// Created by Alex on 19.05.2021.
//

#ifndef IKARUS_OCCUPATIONPATTERN_H
#define IKARUS_OCCUPATIONPATTERN_H

#include <memory>


#include <Eigen/Sparse>

template<typename GridType, typename MatrixType>
class SimpleOccupationPattern
{
//    void generate()
//    void get()


private:
    std::vector< Eigen::Triplet<double> > triplet;
    std::shared_ptr<GridType> grid;
};

#endif //IKARUS_OCCUPATIONPATTERN_H
