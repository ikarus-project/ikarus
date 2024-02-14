// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/parametertree.hh>

#include <ikarus/finiteelements/fesettings.hh>
#include <ikarus/utils/refl.hpp>

namespace Ikarus {

namespace Concepts {
  template <typename T>
  concept Reflectable = refl::is_reflectable<T>();

  template <typename T>
  concept SettingsConatiner = Reflectable<T> and std::is_default_constructible_v<T>;
} // namespace Concepts

/**
 * \brief Settings Template. Can be used with reflectable Settings Containers
 * \tparam C Type of the settings container
 */
template <Concepts::SettingsConatiner C>
struct Settings
{
  using Container = C;

  explicit Settings(const Dune::ParameterTree& parameters)
      : parameters_{parameters} {
    loadSettings();
  }

  /**
   * \brief returns a struct containing the settings
   * \return settings struct
   */
  auto getContainer() const -> const Container& { return settings_; }

  /**
   * \brief Returns the value of a given settings parameter by compile-time string
   * \details
   * \code
   * constexpr auto key = refl::make_const_string("nGP");
   * settings.get<key>();
   * \endcode
   * \tparam str compile time string of parameter name (NB: Will not compile when str is not maching any member)
   * \return value of settings parameter
   */
  template <refl::const_string str>
  constexpr auto& get() const {
    constexpr auto func = find_one(containerReflection_.members, [](auto m) { return m.name == str; });
    return func(settings_);
  }

  /**
   * \brief Helper function to get an easy way to print out settings parameters and values
   */
  friend std::ostream& operator<<(std::ostream& os, const Settings& settings) {
    os << "Settings: " << containerReflection_.name << ": \n";
    for_each(containerReflection_.members, [&]<typename T>(T member) {
      os << "\t" << member.name << ": " << member.get(settings.getContainer()) << " ("
         << Dune::className<typename T::value_type>() << ")\n";
    });
    return os;
  }

private:
  static constexpr auto containerReflection_ = refl::reflect<Container>();
  const Dune::ParameterTree& parameters_;
  Container settings_{};

  auto loadSettings() {
    for_each(containerReflection_.members, [&](auto member) {
      if (member.is_writable)
        parseParameter(static_cast<std::string>(member.name), member.get(settings_));
    });
  }

  template <typename T>
  auto parseParameter(const std::string& key, T& parameter) {
    parameter = parameters_.get<T>(key, parameter);
  }
};

/**
 * Predefined type aliases for Settings
 */
using FESettings = Settings<FESettingsContainer>;

} // namespace Ikarus
