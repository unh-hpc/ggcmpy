
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <openggcm/emfields.hxx>

#include <nanobind/stl/detail/nb_array.h>
#include <xtensor/containers/xfixed.hpp>

NAMESPACE_BEGIN(NB_NAMESPACE)
NAMESPACE_BEGIN(detail)

template <typename Type, size_t Size>
struct type_caster<xt::xtensor_fixed<Type, xt::xshape<Size>>>
    : array_caster<xt::xtensor_fixed<Type, xt::xshape<Size>>, Type, Size>
{
};

NAMESPACE_END(detail)
NAMESPACE_END(NB_NAMESPACE)

namespace nb = nanobind;
using namespace nb::literals;
using namespace openggcm;

NB_MODULE(_openggcm, m)
{
  m.def("xt_fixed_from_python",
        [](xt::xtensor_fixed<double, xt::xshape<3>> a) { return a; });

  nb::class_<emfields>(m, "emfields")
      .def("E", &emfields::E, "r"_a)
      .def("B", &emfields::B, "r"_a)
      .def("__repr__", &emfields::repr);

  nb::class_<emfields_uniform, emfields>(m, "emfields_uniform")
      .def(nb::init<double3, double3>(), "E_0"_a, "B_0"_a);
}
