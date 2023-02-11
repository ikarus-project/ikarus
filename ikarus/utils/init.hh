// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include <dune/common/parallel/mpihelper.hh>
#include <utility>

#include <spdlog/fmt/ostr.h>
#include <spdlog/fmt/ranges.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <spdlog/details/registry.h>

namespace Ikarus {

  class IkarusInstance {
  public:
    static IkarusInstance& getInstance() {
      static IkarusInstance instance;
      return instance;
    }

    void enableFileLogger(std::string&& filename= "") {
      auto logFilename = (filename.empty() ? executableName : filename) + ".log";
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(executableName, true);

        spdlog::default_logger()->sinks().push_back(file_sink);
    }

  private:
    friend void init(int argc, char** argv);
    void setExecutableName(std::string&& p_executableName) { executableName = p_executableName; }
    IkarusInstance() = default;
    std::string executableName;

  public:
    IkarusInstance(IkarusInstance const&) = delete;
    void operator=(IkarusInstance const&) = delete;
  };

  void init(int argc, char** argv) {
    Dune::MPIHelper::instance(argc, argv);
    auto& instance = IkarusInstance::getInstance();
    instance.executableName=argv[0];
  }

  void enableFileLogger(std::string filename= "")
  {
    auto& instance = IkarusInstance::getInstance();
    instance.enableFileLogger(std::move(filename));
  }
}  // namespace Ikarus
