//
// Created by Alex on 21.06.2021.
//

#pragma once

#include <map>

#include <ikarus/FiniteElements/FiniteElementInterface.h>
#include <ikarus/FiniteElements/GenericFiniteElement.h>
/*
 * See API design for C++, Martin Reddy, Chapter 3.3
 */
namespace Ikarus::FiniteElements {

  class FiniteElementFactory {
  public:
    template <typename FEType> static void registerFiniteElement(const std::string& type, FEType fe);
    static void unRegisterFiniteElement(const std::string& type);

  private:
    using FEMap = std::map<std::string, GenericFE>;
    static FEMap femap;
  };
}  // namespace Ikarus::FiniteElements