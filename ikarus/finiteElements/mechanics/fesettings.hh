// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <map>
#include <string>

namespace Ikarus
{
  struct FESettings{

    template<typename T>
    void addOrAssign(const std::string& key , T value)
    {
      if constexpr (std::is_same_v<T,bool>)
      boolMap.insert_or_assign(key,value);
      else  if constexpr (std::is_same_v<T,double>)
      doubleMap.insert_or_assign(key,value);
            else  if constexpr (std::is_same_v<T,int>)
      intMap.insert_or_assign(key,value);
    }

    template<typename T>
    auto request(const std::string& key ) const 
    {
      if constexpr (std::is_same_v<T,bool>)
      return boolMap.at(key);
      else  if constexpr (std::is_same_v<T,double>)
      return doubleMap.at(key);
             else  if constexpr (std::is_same_v<T,int>)
      return intMap.at(key);
    }


    std::map<std::string,bool> boolMap;
    std::map<std::string,double> doubleMap;
    std::map<std::string,int > intMap;


  };
} // namespace Ikarus
