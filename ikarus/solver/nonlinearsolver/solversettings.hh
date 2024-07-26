// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file solversettings.hh
 * \brief Macros to define solver setting structs
 */

#pragma once



#define MEMBERVARIABLES(name, type, default_value, comment) \
    type name{default_value}; /* #comment */

#define MEMBERVARIABLESSETTER(name, type, default_value, comment)       \
    {#name, [this](const pybind11::handle &value) {                           \
         try {                                                          \
             name = value.cast<type>();                                 \
                                                                        \
         } catch (const pybind11::cast_error &e) {                            \
             pybind11::cast_error eAmmend(                                    \
                 ("The value type of \"" + std::string(#name) +         \
                  "\" should be something as <" + std::string(#type) +  \
                  "> "                                                  \
                  "but the cast failed. \nThe original Pybind11 error " \
                  "reads:\n") +                                         \
                 e.what());                                             \
             throw eAmmend;                                             \
         }                                                              \
     }},

#define SOLVERSETTINGS(Name, FIELDS)                                  \
    struct Name {                                                     \
        FIELDS(MEMBERVARIABLES)                                       \
                                                                      \
        void populate(const pybind11::dict &dict) {                         \
            static const std::unordered_map<                          \
                std::string, std::function<void(const pybind11::handle &)>> \
                setters = {FIELDS(MEMBERVARIABLESSETTER)};\
\
for (auto item : dict) {\
    std::string key = item.first.cast<std::string>();\
    if (setters.count(key) > 0) {\
        setters.at(key)(item.second);\
    } else\
        throw std::runtime_error("Invalid key in dictionary: " + key);\
}}};
