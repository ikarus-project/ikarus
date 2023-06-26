// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

namespace Ikarus
{
    template <typename Element, typename NeumannBoundary, typename F>
    auto forEachInterSectionIntegrationPoint(const Element& element, const NeumannBoundary& neumannBoundary, int order, F&& f)
    {
      for (auto&& intersection : intersections(neumannBoundary->gridView(), element)) {
        if (not neumannBoundary or not neumannBoundary->contains(intersection)) continue;

        const auto& quadLine = Dune::QuadratureRules<double, Element::mydimension - 1>::rule(intersection.type(), order);

        for (const auto& curQuad : quadLine) {
          // Local position of the quadrature point
          const Dune::FieldVector<double,  Element::mydimension>& quadPos
              = intersection.geometryInInside().global(curQuad.position());

          const double intElement = intersection.geometry().integrationElement(curQuad.position())*curQuad.weight();
          f(quadPos,Dune::toEigen(intersection.geometry().global(curQuad.position())),intElement) ;
        }
      }
    }


}
