//
// Created by alex on 2/8/22.
//

#pragma once
#include <ranges>
 namespace Ikarus
{
   template<typename Basis, typename Predicate >
   void markDirichletBoundaryDofs(const Basis& basis,std::ranges::random_access_range  auto& dirichletFlags, Predicate&& f)
   {
     auto seDOFs = subEntityDOFs(basis);
     const auto& gridView = basis.gridView();
     auto localView = basis.localView();
     for (auto&& element : elements(basis.gridView()))
       if (element.hasBoundaryIntersections()) {
         localView.bind(element);
         for (const auto& intersection : intersections(gridView, element))
           if (intersection.boundary())
             if (f(intersection.geometry().center()))
               for (auto localIndex : seDOFs.bind(localView, intersection))
                 dirichletFlags[localView.index(localIndex)[0]] = true;
       }
   }
 }