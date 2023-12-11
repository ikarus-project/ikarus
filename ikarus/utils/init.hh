// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include <chrono>
#include <utility>

#include <dune/common/parallel/mpihelper.hh>

#include <spdlog/details/registry.h>
#include <spdlog/fmt/chrono.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/fmt/ranges.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace Ikarus {

  class IkarusInstance {
  public:
    static IkarusInstance& getInstance() {
      static IkarusInstance instance;
      return instance;
    }

    void enableFileLogger(std::string&& filename = "") {
      using namespace std::chrono;
      std::string currentTime = fmt::format("_{}", std::chrono::system_clock::now());

      std::ranges::transform(currentTime, currentTime.begin(), [](char ch) {
        return (ch == ' ' or ch == ':') ? '_' : ch;
      });  // replace space and colon with underscore
      auto logFilename = (filename.empty() ? executableName : filename) + currentTime + ".log";
      auto file_sink   = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logFilename, true);
      file_sink->set_pattern("%v");
      file_sink->set_level(spdlog::level::trace);
      auto logger = spdlog::default_logger();
      logger->sinks().push_back(file_sink);
    }

  private:
    friend void init(int argc, char** argv, bool enableFileLogger);
    IkarusInstance() = default;
    std::string executableName;

  public:
    IkarusInstance(IkarusInstance const&) = delete;
    void operator=(IkarusInstance const&) = delete;
  };

  void init(int argc, char** argv, bool enableFileLogger = true) {
    Dune::MPIHelper::instance(argc, argv);
    auto& instance          = IkarusInstance::getInstance();
    instance.executableName = argv[0];
    auto logger             = spdlog::default_logger();
    logger->set_pattern("%v");
    if (enableFileLogger) instance.enableFileLogger();

    auto currentTime = std::chrono::system_clock::now();
    std::chrono::year_month_day currentYearMonthDay{floor<std::chrono::days>(currentTime)};

    spdlog::info("Start of execution: {}", currentTime);
    /// https://patorjk.com/software/taag/#p=testall&f=Univers&t=IKARUS (font: Lean)
    spdlog::info(R"xxx(    _/_/_/  _/    _/    _/_/    _/_/_/    _/    _/    _/_/_/)xxx");
    spdlog::info(R"xxx(     _/    _/  _/    _/    _/  _/    _/  _/    _/  _/       )xxx");
    spdlog::info(R"xxx(    _/    _/_/      _/_/_/_/  _/_/_/    _/    _/    _/_/    )xxx");
    spdlog::info(R"xxx(   _/    _/  _/    _/    _/  _/    _/  _/    _/        _/   )xxx");
    spdlog::info(R"xxx(_/_/_/  _/    _/  _/    _/  _/    _/    _/_/    _/_/_/      )xxx");

    spdlog::info("© 2021-{} The Ikarus Developers, see LICENSE.md ", static_cast<int>(currentYearMonthDay.year()));
    spdlog::info("You are using Ikarus v{}. Please don't forget to cite us:", IKARUS_VERSION);
    spdlog::info(
        "Müller, A., & Vinod Kumar Mitruka, T. K. M. (2023). Ikarus v0.3 (Version V1). Version V1. "
        "doi:<https://doi.org/10.18419/darus-3303>");
  }
}  // namespace Ikarus
