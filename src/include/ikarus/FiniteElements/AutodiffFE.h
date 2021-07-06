////
//// Created by Alex on 05.07.2021.
////
//
//#pragma once
////This element can not be used on its own but it should be inherited from
//// The class constructor can only be called from the templated class.
//template<typename ConcreteElement>
//class AutoDiffFE
//{
//  constexpr Derived&       derived()       noexcept { return static_cast<Derived&>( *this ); }
//  constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>( *this ); }
//
//  constexpr adouble adEnergy() const noexcept
//  {
//    return derived().energyImpl<adouble>(); }
//
//private:
//  AutoDiffFE() = default;
//  friend class ConcreteElement;
//};