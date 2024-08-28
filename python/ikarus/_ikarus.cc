// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "io/io.hh"
#include "pythonhelpers.hh"

#include <dune/common/float_cmp.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <spdlog/spdlog.h>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/python/finiteelements/material.hh>
#include <ikarus/python/finiteelements/scalarwrapper.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>

/**
 * \brief Registers the ScalarWrapper class template with pybind11, adding various operations and constructors.
 *
 * This function template registers the ScalarWrapper class in the given Python module with specified names for
 * the class and its scalar type. It adds support for basic arithmetic operations, comparisons, and string
 * representation.
 *
 * \tparam ScalarType The type of the scalar value to be wrapped.
 * \param m The pybind11 module where the class will be registered.
 * \param name The name of the class in the Python module.
 * \param typeName The cpp name of the scalar type in the Python module.
 */
template <typename ScalarType>
void registerScalarWrapper(pybind11::module& m, std::string name, std::string typeName) {
  namespace py = pybind11;
  using namespace Ikarus;
  using namespace Eigen;
  using VWd           = ScalarWrapper<ScalarType>;
  using RawScalarType = Dune::ResolveRef_t<ScalarType>;

  auto scalarWrapper = Dune::Python::insertClass<VWd>(m, name, Dune::Python::GenerateTypeName(typeName)).first;
  if constexpr (std::is_same_v<VWd, RawScalarType>)
    scalarWrapper.def(py::init<RawScalarType>(), "Constructor with a raw scalar value.");
  scalarWrapper.def(py::init<VWd>(), "Copy constructor.");
  scalarWrapper.def(py::init<ScalarType>(), "Constructor with a scalar value.");
  scalarWrapper.def("__repr__", [](const VWd& d) { return std::to_string(d.value()); }, "String representation.");
  scalarWrapper.def(
      "__eq__", [](const VWd& x, const VWd& y) { return Dune::FloatCmp::eq(x.value(), y.value()); },
      "Equality comparison with another ScalarWrapper.");
  scalarWrapper.def(
      "__eq__", [](const VWd& x, const RawScalarType& y) { return Dune::FloatCmp::eq(x.value(), y); },
      "Equality comparison with a raw scalar value.");
  scalarWrapper.def(
      "__eq__", [](const RawScalarType& x, const VWd& y) { return Dune::FloatCmp::eq(x, y.value()); },
      "Equality comparison of a raw scalar value with a ScalarWrapper.");
  scalarWrapper.def(py::self + py::self, "Addition of two ScalarWrapper instances.")
      .def(py::self - py::self, "Subtraction of two ScalarWrapper instances.")
      .def(py::self += py::self, "Addition assignment of another ScalarWrapper instance.")
      .def(py::self -= py::self, "Subtraction assignment of another ScalarWrapper instance.")
      .def(-py::self, "Negation of the ScalarWrapper value.")
      .def(py::self *= RawScalarType(), "Multiplication assignment with a raw scalar value.")
      .def(py::self /= RawScalarType(), "Division assignment with a raw scalar value.")
      .def(RawScalarType() * py::self, "Multiplication of a raw scalar value with a ScalarWrapper.")
      .def(py::self * RawScalarType(), "Multiplication of a ScalarWrapper with a raw scalar value.")
      .doc() =
      R"(Since Python does not support passing python float by reference to a double&, we have to wrap everything to allow mutable scalar types that can be passed back and forth.
              See also: https://pybind11.readthedocs.io/en/stable/faq.html#limitations-involving-reference-arguments)";
}

PYBIND11_MODULE(_ikarus, m) {
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;
  using namespace Eigen;

  ENUM_BINDINGS(ScalarAffordance);
  ENUM_BINDINGS(VectorAffordance);
  ENUM_BINDINGS(MatrixAffordance);
  ENUM_BINDINGS(FESolutions);
  ENUM_BINDINGS(FEParameter);
  ENUM_BINDINGS(SolverTypeTag);
  ENUM_BINDINGS(DBCOption);

  py::class_<std::tuple<ScalarAffordance, VectorAffordance, MatrixAffordance>>(m, "Base");
  auto affordanceCollections3 =
      Dune::Python::insertClass<AffordanceCollection<ScalarAffordance, VectorAffordance, MatrixAffordance>,
                                std::tuple<ScalarAffordance, VectorAffordance, MatrixAffordance>>(
          m, "AffordanceCollection",
          Dune::Python::GenerateTypeName("AffordanceCollection<ScalarAffordance,VectorAffordance,MatrixAffordance>"),
          Dune::Python::IncludeFiles("ikarus/finiteelements/ferequirements.hh"))
          .first;

  affordanceCollections3.def(py::init<ScalarAffordance, VectorAffordance, MatrixAffordance>());
  affordanceCollections3.def_property_readonly_static("elastoStatics",
                                                      [](py::object) { return AffordanceCollections::elastoStatics; });

  registerScalarWrapper<double>(m, "Scalar", "ScalarWrapper<double>");
  registerScalarWrapper<std::reference_wrapper<double>>(m, "ScalarRef",
                                                        "ScalarWrapper<std::reference_wrapper<double>>");

  // This should maybe be moved to a different .cc file for the submodule
  auto materials = m.def_submodule("materials", "This is the submodule for materials in Ikarus");

  ENUM_BINDINGS_WITH_MODULE(StrainTags, materials);
  ENUM_BINDINGS_WITH_MODULE(StressTags, materials);
  ENUM_BINDINGS_WITH_MODULE(TangentModuliTags, materials);

  pybind11::class_<LinearElasticity> linElastic(materials, "LinearElasticity");
  Ikarus::Python::registerLinearElasticity(materials, linElastic);

  pybind11::class_<StVenantKirchhoff> svk(materials, "StVenantKirchhoff");
  Ikarus::Python::registerStVenantKirchhoff(materials, svk);

  pybind11::class_<NeoHooke> nh(materials, "NeoHooke");
  Ikarus::Python::registerNeoHooke(materials, nh);

  // This could go into a submodule util
  materials.def(
      "toVoigt",
      [](Eigen::MatrixXd mat, bool isStrain = true) {
        auto callToVoigt = []<typename T>(const T& m, bool isStrain) {
          auto result = toVoigt(m, isStrain);
          return Eigen::MatrixXd(result);
        };

        if (mat.rows() == 1 && mat.cols() == 1) {
          Eigen::Matrix<double, 1, 1> fixedMat = mat;
          return callToVoigt(fixedMat, isStrain);
        } else if (mat.rows() == 2 && mat.cols() == 2) {
          Eigen::Matrix<double, 2, 2> fixedMat = mat;
          return callToVoigt(fixedMat, isStrain);

        } else if (mat.rows() == 3 && mat.cols() == 3) {
          Eigen::Matrix<double, 3, 3> fixedMat = mat;
          return callToVoigt(fixedMat, isStrain);
        } else {
          DUNE_THROW(Dune::IOError, "toVoigt only supports matrices of dimension 1, 2 or 3");
        }
      },
      py::arg("matrix"), py::arg("isStrain") = true);

  materials.def(
      "fromVoigt",
      [](Eigen::MatrixXd vec, bool isStrain = true) {
        auto callFromVoigt = []<typename T>(const T& m, bool isStrain) {
          auto result = fromVoigt(m, isStrain);
          return Eigen::MatrixXd(result);
        };
        if (vec.rows() == 1) {
          Eigen::Vector<double, 1> fixedMat = vec;
          return callFromVoigt(fixedMat, isStrain);
        } else if (vec.rows() == 3) {
          Eigen::Vector<double, 3> fixedMat = vec;
          return callFromVoigt(fixedMat, isStrain);

        } else if (vec.rows() == 6) {
          Eigen::Vector<double, 6> fixedMat = vec;
          return callFromVoigt(fixedMat, isStrain);
        } else {
          DUNE_THROW(Dune::IOError, "fromVoigt only supports matrices of dimension 1, 2 or 3");
        }
      },
      py::arg("vector"), py::arg("isStrain") = true);

  /**
   * \brief Transform strain from one type to another.
   *
   * This function transforms one strain component matrix from one type to another, based on the provided strain tags
   * \ingroup  materials
   * \param from Type of the source strain tag.
   * \param to Type of the target strain tag.
   * \param E Eigen matrix representing the input strain (can be in Voigt notation).
   * \return The transformed strain matrix.
   */
  materials.def(
      "tramsformStrain",
      [](StrainTags from, StrainTags to, Eigen::MatrixXd E) {
        auto callTransformStrain =
            []<StrainTags from_, StrainTags to_>(const Eigen::MatrixXd& eRaw) -> Eigen::MatrixXd {
          return transformStrain<from_, to_>(Eigen::Matrix<double, 3, 3>(eRaw));
        };
        if ((E.rows() != 3 and E.cols() != 3) or (E.rows() == 6))
          DUNE_THROW(Dune::IOError,
                     "Strain converseions are only implemented for matrices of dimension 3 or the corresponding Voigt "
                     "notation");

        if (from == StrainTags::linear) {
          spdlog::warn("No useful transformation available for linear strains");
          return E;
        }
        if (from == to)
          return E;

        if (to == StrainTags::greenLagrangian) {
          switch (from) {
            case StrainTags::deformationGradient:
              return callTransformStrain
                  .template operator()<StrainTags::deformationGradient, StrainTags::greenLagrangian>(E);
            case Ikarus::StrainTags::displacementGradient:
              return callTransformStrain
                  .template operator()<StrainTags::displacementGradient, StrainTags::greenLagrangian>(E);
            case Ikarus::StrainTags::rightCauchyGreenTensor:
              return callTransformStrain
                  .template operator()<StrainTags::rightCauchyGreenTensor, StrainTags::greenLagrangian>(E);
            default:
              __builtin_unreachable();
          }
        } else if (to == StrainTags::deformationGradient) {
          switch (from) {
            case StrainTags::greenLagrangian:
              return callTransformStrain
                  .template operator()<StrainTags::greenLagrangian, StrainTags::deformationGradient>(E);
            case Ikarus::StrainTags::displacementGradient:
              return callTransformStrain
                  .template operator()<StrainTags::displacementGradient, StrainTags::deformationGradient>(E);
            case Ikarus::StrainTags::rightCauchyGreenTensor:
              return callTransformStrain
                  .template operator()<StrainTags::rightCauchyGreenTensor, StrainTags::deformationGradient>(E);
            default:
              __builtin_unreachable();
          }
        } else if (to == StrainTags::rightCauchyGreenTensor) {
          switch (from) {
            case StrainTags::greenLagrangian:
              return callTransformStrain
                  .template operator()<StrainTags::greenLagrangian, StrainTags::rightCauchyGreenTensor>(E);
            case Ikarus::StrainTags::displacementGradient:
              return callTransformStrain
                  .template operator()<StrainTags::displacementGradient, StrainTags::rightCauchyGreenTensor>(E);
            case Ikarus::StrainTags::deformationGradient:
              return callTransformStrain
                  .template operator()<StrainTags::deformationGradient, StrainTags::rightCauchyGreenTensor>(E);
            default:
              __builtin_unreachable();
          }
        }
        __builtin_unreachable();
      },
      py::arg("to"), py::arg("from"), py::arg("strain"));

  addBindingsToIO();
}
