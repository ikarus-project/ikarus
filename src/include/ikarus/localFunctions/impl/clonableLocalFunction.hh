//
// Created by lex on 5/1/22.
//

#pragma once
#include <ikarus/localFunctions/meta.hh>
namespace Ikarus {

  template<typename LFImpl>
  class ClonableLocalFunction {

  public:
    LFImpl clone() const {
        return LFImpl(underlying().basis(),underlying().coeffs,typename LFImpl::Ids() );
    }

    template <typename OtherType,size_t ID=0>
    auto rebindClone(OtherType&& t,Dune::index_constant<ID> && id = Dune::index_constant<0>() ) const
    {
      if constexpr (LFImpl::Ids::value == ID)
        return typename LFImpl::template Rebind<OtherType>::other(underlying().basis(), convertUnderlying<OtherType>(underlying().coeffs),typename LFImpl::Ids() );
      else
        return clone();
    }



  private:
    constexpr LFImpl const& underlying() const  // CRTP
    {
      return static_cast<LFImpl const&>(*this);
    }


  };

}  // namespace Ikarus

