
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>

#include <openggcm/emfields.hxx>

#include <xtensor/containers/xfixed.hpp>
#include <xtensor/containers/xadapt.hpp>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_openggcm, m)
{
  m.def("xt_fixed_from_python", [](std::array<double, 3> a) {
    auto r = xt::xtensor_fixed<double, xt::xshape<3>>{a[0], a[1], a[2]};
    // auto r = xt::adapt(a.data(), {a.shape(0)}); // xt::xshape<3>());
    return nb::ndarray<nb::numpy, const double, nb::shape<3>>(r.data()).cast();
  });

  nb::class_<openggcm::emfields_uniform>(m, "emfields_uniform")
    .def(nb::init<std::array<double, 3>, std::array<double, 3>>(), "E_0"_a,
         "B_0"_a)
    .def(
      "E",
      [](openggcm::emfields_uniform& self, std::array<double, 3> r) {
        return self.E(r[0], r[1], r[2]);
      },
      "r"_a)
    .def(
      "B",
      [](openggcm::emfields_uniform& self, std::array<double, 3> r) {
        return self.B(r[0], r[1], r[2]);
      },
      "r"_a)
    .def("__repr__", &openggcm::emfields_uniform::repr);
}
