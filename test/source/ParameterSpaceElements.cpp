//
// Created by Alex on 25.05.2021.
//
#include "gtest/gtest.h"

#include <dune/geometry/referenceelements.hh>

#include <ikarus/ParameterSpaceElements/LinearCube.h>

TEST(ParaSpaceElementTest, Test) {
  Ikarus::ElementTypes::LinearCube<1> c1;
  Ikarus::ElementTypes::LinearCube<2> c2;
  Ikarus::ElementTypes::LinearCube<3> c3;

  EXPECT_EQ(c2.numberOfEdges(), 4);
  EXPECT_EQ(c3.numberOfEdges(), 12);

  EXPECT_EQ(c1.dim(), 1);
  EXPECT_EQ(c2.dim(), 2);
  EXPECT_EQ(c3.dim(), 3);

  EXPECT_EQ(c1.numberOfVertices(), 2);
  EXPECT_EQ(c2.numberOfVertices(), 4);
  EXPECT_EQ(c3.numberOfVertices(), 8);

  EXPECT_EQ(c3.numberOfSurfaces(), 6);
}