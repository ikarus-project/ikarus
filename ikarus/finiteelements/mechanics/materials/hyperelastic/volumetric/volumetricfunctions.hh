// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file volumetricfunctions.hh
 * \brief Implementation of the volumetric part of a hyperelastic material.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/finiteelements/mechanics/materials/hyperelastic/volumetric/interface.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Default volumetric function
 * \details \f[ U(J) = 0 \f]
 * \ingroup materials
 */

struct VF0
{
  template <typename ST>
  ST storedEnergyImpl(const ST& /* J */) const {
    return 0;
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& /* J */) const {
    return 0;
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& /* J */) const {
    return 0;
  };

  [[nodiscard]] constexpr static std::string name() noexcept { return "None"; }
};

/**
 * \brief Volumetric function No. 1 found in \cite hartmann_polyconvexity_2003 Tab. 4
 * \details \f[ U(J) = \frac{1}{2}(J - 1)^2 \f].
 * \ingroup materials
 */

struct VF1
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return 0.5 * pow(J - 1, 2);
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return J - 1;
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& /* J */) const {
    return 1;
  };

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 1"; }
};

/**
 * \brief Volumetric function No. 2 found in \cite hartmann_polyconvexity_2003 Tab. 4
 * \details \f[ U(J) = \frac{1}{4}((J - 1)^2 + (\ln J )^2) \f].
 * \ingroup materials
 */

struct VF2
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return 0.25 * (pow(J - 1, 2) + pow(log(J), 2));
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return 0.5 * (J - 1 + 1 / J * log(J));
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    auto Jsq = pow(J, 2);
    return 1 / (2 * Jsq) * (1 + Jsq - log(J));
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 2"; }
};

/**
 * \brief Volumetric function No. 3 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = \frac{1}{2}(\ln J )^2 \f]
 * \ingroup materials
 */

struct VF3
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return 0.5 * pow(log(J), 2);
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return 1 / J * log(J);
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return 1 / pow(J, 2) * (1 - log(J));
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 3"; }
};

/**
 * \brief Volumetric function No. 4 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = \frac{1}{\beta^2}(\frac{1}{J^\beta} - 1 + \beta(\ln J)) \f]
 * \ingroup materials
 */

struct VF4
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return (1 / pow(beta_, 2)) * ((1 / pow(J, beta_)) - 1 + beta_ * log(J));
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return (1 / beta_) * ((1 / J) - (1 / (pow(J, 1 + beta_))));
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return (1 / beta_) * ((1 / pow(J, 2 + beta_)) * (1 + beta_ - pow(J, beta_)));
  }

  explicit VF4(double beta)
      : beta_(beta) {}

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 4"; }

  double beta() const { return beta_; }

private:
  double beta_;
};

/**
 * \brief Volumetric function No. 5 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = \frac{1}{4}(J^2 - 1 - 2 \ln J) \f]
 * \ingroup materials
 */

struct VF5
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return 0.25 * (pow(J, 2) - 1 - 2 * log(J));
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return 0.5 * (J - (1 / J));
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return 0.5 * (1 + (1 / pow(J, 2)));
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 5"; }
};

/**
 * \brief Volumetric function No. 6 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = (J - \ln J - 1) \f]
 * \ingroup materials
 */

struct VF6
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return J - log(J) - 1;
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return 1 - (1 / J);
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return 1 / pow(J, 2);
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 6"; }
};

/**
 * \brief Volumetric function No. 7 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = J^\beta(\beta \ln J - 1) + 1 \f]
 * \ingroup materials
 */

struct VF7
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return pow(J, beta_) * (beta_ * log(J) - 1) + 1;
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return pow(beta_, 2) * (1.0 / pow(J, 1.0 - beta_)) * log(J);
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return pow(beta_, 2) * pow(J, beta_ - 2.0) * (1 + (beta_ - 1) * log(J));
  }

  explicit VF7(double beta)
      : beta_(beta) {}

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 7"; }

  double beta() const { return beta_; }

private:
  double beta_;
};

/**
 * \brief Volumetric function No. 8 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = J \ln J - J + 1 \f].
 * \ingroup materials
 */

struct VF8
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return J * log(J) - J + 1;
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return log(J);
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return 1 / J;
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 8"; }
};

/**
 * \brief Volumetric function No. 9 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = \frac{1}{32}(J^2 - J^{-2})^2 \f].
 * \ingroup materials
 */

struct VF9
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return (1.0 / 32.0) * pow(pow(J, 2) - pow(J, -2), 2);
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return (1.0 / 8.0) * (pow(J, 3) - (1.0 / pow(J, 5)));
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return (1.0 / 8.0) * (5.0 * pow(J, -6) + (3.0 * pow(J, 2)));
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 9"; }
};

/**
 * \brief Volumetric function No. 10 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = \frac{J}{\beta}(1 - \frac{J^{-\beta}}{1-\beta}) +
 \frac{1}{\beta - 1} \f]
 * \ingroup materials
 */

struct VF10
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return (J / beta_) * (1 - (pow(J, -beta_) / (1 - beta_))) + (1.0 / (beta_ - 1));
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return (1 / beta_) * (1 - pow(J, -beta_));
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return pow(J, -1 - beta_);
  }

  explicit VF10(double beta)
      : beta_(beta) {}

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 10"; }

  double beta() const { return beta_; }

private:
  double beta_;
};

/**
 * \brief Volumetric function No. 11 found in \cite hartmann_polyconvexity_2003 Tab. 4.
 * \details \f[ U(J) = \frac{1}{50}(J^5 + J^{-5} - 2) \f]
 * \ingroup materials
 */

struct VF11
{
  template <typename ST>
  ST storedEnergyImpl(const ST& J) const {
    return (1.0 / 50.0) * (pow(J, 5.0) + pow(J, -5.0) - 2.0);
  };

  template <typename ST>
  ST firstDerivativeImpl(const ST& J) const {
    return (1.0 / 10.0) * (pow(J, 4.0) - pow(J, -6.0));
  }

  template <typename ST>
  ST secondDerivativeImpl(const ST& J) const {
    return (1.0 / 10.0) * (4 * pow(J, 3.0) + 6 * pow(J, -7.0));
  }

  [[nodiscard]] constexpr static std::string name() noexcept { return "Function 11"; }
};

using NoVolumetricPart = Volumetric<VF0>;

} // namespace Ikarus::Materials
