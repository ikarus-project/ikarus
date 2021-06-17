//
// Created by Alex on 27.05.2021.
//

#pragma once

/** \brief A traits which returns false if template is instiantiated, handy für static_assert */
template<typename... Args>
bool dependentFalse()
{
    return false;
}