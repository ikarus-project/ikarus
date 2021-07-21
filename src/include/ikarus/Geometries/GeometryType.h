#pragma once
#include <iostream>

// Further types should be added in the future
namespace Ikarus{
  enum class GeometryType {
    vertex,
    linearLine,
    quadraticLine,
    cubicLine,
    linearTriangle,
    linearQuadrilateral,
    linearTetrahedron,
    pyramid,
    prism,
    linearHexahedron
  };

  inline bool isVertex(GeometryType type){return type == Ikarus::GeometryType::vertex;}
  inline bool isLinearLine(GeometryType type){return type == Ikarus::GeometryType::linearLine;}
  inline bool isQuadraticLine(GeometryType type){return type == Ikarus::GeometryType::quadraticLine;}
  inline bool isCubicLine(GeometryType type){return type == Ikarus::GeometryType::cubicLine;}
  inline bool isLinearTriangle(GeometryType type){return type == Ikarus::GeometryType::linearTriangle;}
  inline bool isLinearQuadrilateral(GeometryType type){return type == Ikarus::GeometryType::linearQuadrilateral;}
  inline bool isLinearTetrahedron(GeometryType type){return type == Ikarus::GeometryType::linearTetrahedron;}
  inline bool isPyramid(GeometryType type){return type == Ikarus::GeometryType::pyramid;}
  inline bool isPrism(GeometryType type){return type == Ikarus::GeometryType::prism;}
  inline bool isLinearHexahedron(GeometryType type){return type == Ikarus::GeometryType::linearHexahedron;}

  inline int dimension(GeometryType type)
  {
    if (isVertex(type)){
      return 0;
    }
    else if (isLinearLine(type) || isQuadraticLine(type) || isCubicLine(type)){
      return 1;
    }
    else if (isLinearTriangle(type) || isLinearQuadrilateral(type)){
      return 2;
    }
    else if (isLinearTetrahedron(type) || isPyramid(type) || isPrism(type) || isLinearHexahedron(type)){
      return 3;
    }
    else {
      throw std::runtime_error("Dimension of the given type couldn't be obtained. Check function 'dimension' in GeometryType.h");
    }
  }



  inline std::string print (GeometryType type)
  {
    switch (type)
    {
      case GeometryType::vertex : return "vertex";
      case GeometryType::linearLine : return "linearLine";
      case GeometryType::quadraticLine : return "quadraticLine";
      case GeometryType::cubicLine : return "cubicLine";
      case GeometryType::linearTriangle : return "linearTriangle";
      case GeometryType::linearQuadrilateral : return "linearQuadrilateral";
      case GeometryType::linearTetrahedron : return "linearTetrahedron";
      case GeometryType::pyramid : return "pyramid";
      case GeometryType::prism : return "prism";
      case GeometryType::linearHexahedron : return "linearHexahedron";
      default: throw std::runtime_error("ostream operator for this GeometryType is not implemented. Check implementation in GeometryType.h");
    }
  }

  inline std::ostream& operator<< (std::ostream& os, GeometryType type)
  {
    os << print(type);
    return os;
  }

  // This function is necessary to work with Ikarus Geometry-Types in the GridEntity implementation based on Dune
  inline Dune::GeometryType duneType(Ikarus::GeometryType ikarusType)
  {
    switch (ikarusType)
    {
      case GeometryType::vertex : return Dune::GeometryTypes::vertex;
      case GeometryType::linearLine : return Dune::GeometryTypes::line;
      case GeometryType::linearTriangle : return Dune::GeometryTypes::triangle;
      case GeometryType::linearQuadrilateral : return Dune::GeometryTypes::quadrilateral;
      case GeometryType::linearTetrahedron : return Dune::GeometryTypes::tetrahedron;
      case GeometryType::pyramid : return Dune::GeometryTypes::pyramid;
      case GeometryType::prism : return Dune::GeometryTypes::prism;
      case GeometryType::linearHexahedron : return Dune::GeometryTypes::hexahedron;
      default: throw std::runtime_error("The type " + print(ikarusType) +" has no equivalent Dune::GeometryType.");
    }
  }

}