//
// Created by Alex on 16.06.2021.
//

#pragma once
#include <ikarus/Grids/GridEntities/GridEntitiesInterface.h>

#include <ikarus/Variables/GenericVariable.h>
#include <ikarus/Variables/VariableInterface.h>

template<Ikarus::Concepts::GridEntity GridEntityType>
class DOFOwnerAdapter {
 public:

  using ctype = GridEntityType::ctype;

  DOFOwnerAdapter(GridEntityType& gE)
  : gridEntity{gE}

  void addDOF(Ikarus::Variable::GenericVariable auto&& var)
  {
    variableVector.push_back(var);
  }




 private:
  std::vector<Ikarus::Variable::GenericVariable> variableVector;
  GridEntityType* gridEntity;
};



