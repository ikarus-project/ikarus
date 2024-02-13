// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/parametertree.hh>

#include <ikarus/finiteelements/fesettings.hh>
#include <ikarus/utils/refl.hpp>

namespace Ikarus {

template <typename FE>
struct Settings
{
  using Container = SettingsContainer;

  explicit Settings(const Dune::ParameterTree& parameters)
      : parameters_{parameters} {
    loadSettings();
  }

  auto getContainer() const -> const Container& { return settings_; }

  template <refl::const_string str>
  constexpr auto& get() const {
    constexpr auto func = find_one(settingsReflection.members, [](auto m) { return m.name == str; });
    return func(settings_);
  }

private:
  static constexpr auto settingsReflection = refl::reflect<Container>();
  const Dune::ParameterTree& parameters_;
  Container settings_{};

  auto loadSettings() {
    for_each(settingsReflection.members, [&](auto member) {
      if (member.is_writable)
        parseParameter(static_cast<std::string>(member.name), member(settings_));
    });
  }

  template <typename T>
  auto parseParameter(const std::string& key, T& parameter) {
    parameter = parameters_.get<T>(key, parameter);
  }
};
} // namespace Ikarus
