
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>

#include <openggcm/emfields.hxx>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_openggcm, m)
{
  nb::class_<openggcm::emfields_constant>(m, "emfields_constant")
    .def(nb::init<std::array<double, 3>, std::array<double, 3>>(), "E0"_a,
         "B0"_a)
    .def("E", &openggcm::emfields_constant::E, "x"_a, "y"_a, "z"_a)
    .def("B", &openggcm::emfields_constant::B, "x"_a, "y"_a, "z"_a)
    .def("__repr__", &openggcm::emfields_constant::repr);
}
