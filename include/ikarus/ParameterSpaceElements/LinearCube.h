//
// Created by Alex on 27.05.2021.
//

#pragma once
#include <dune/geometry/referenceelements.hh>

namespace Ikarus{

    namespace Impl{

    }


namespace ElementTypes
{
    template<int mydim>
    class LinearCube
    {
    public:
        static constexpr int dim(){return duneReferenceElement.dimension;}

//        template<int dummydim = mydim, std::enable_if_t<(dummydim>2) ,bool> = true >
        static constexpr int numberOfSurfaces(){
            static_assert(mydim>2,"Cube dimension must be larger than 2 to have surfaces as subentities.");
            return duneReferenceElement.size(mydim-2);
        }

//        template<int dummydim = mydim, std::enable_if_t<(dummydim>1),bool> = true >
        static constexpr int numberOfEdges(){
            static_assert(mydim>1,"Cube dimension must be larger than 1 to have edges as subentities.");
            return duneReferenceElement.size(mydim-1);
        }
        static constexpr int numberOfVertices(){return duneReferenceElement.size(mydim);}

    private:
        inline static const auto duneReferenceElement = Dune::ReferenceElements<double,mydim>::cube();
    };

    inline constexpr auto q1_1d = LinearCube<1>();
    inline constexpr auto q1_2d = LinearCube<2>();
    inline constexpr auto q1_3d = LinearCube<3>();
//    inline constexpr auto q2 = QuadraticCube<2>();
//    inline constexpr auto p1 = LinearPrism<2>();
//    inline constexpr auto p2 = QuadraticPrism<2>();
}
}

