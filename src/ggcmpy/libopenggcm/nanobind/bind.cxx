
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

namespace test
{
class test_ndarray
{
public:
  test_ndarray(nb::ndarray<double, nb::ndim<1>> arr) : arr_(arr) {}

  double operator[](std::size_t i) const { return arr_(i); }

private:
  nb::ndarray<double, nb::ndim<1>> arr_;
};

class test_xtadapt
{
public:
  test_xtadapt(nb::ndarray<double, nb::ndim<1>> a)
      : a_(std::make_unique<nb::ndarray<double, nb::ndim<1>>>(a))
  {
  }

  double operator[](std::size_t i) const { return a_.get()->operator()(i); }

private:
  std::unique_ptr<nb::ndarray<double, nb::ndim<1>>> a_;
};
} // namespace test

NB_MODULE(_openggcm, m)
{
  m.def("xt_fixed_from_python",
        [](xt::xtensor_fixed<double, xt::xshape<3>> a) { return a; });

  nb::class_<test::test_ndarray>(m, "test_ndarray")
      .def(nb::init<nb::ndarray<double, nb::ndim<1>>>())
      .def("__getitem__", &test::test_ndarray::operator[], "i"_a);

  nb::class_<test::test_xtadapt>(m, "test_xtadapt")
      .def(nb::init<nb::ndarray<double, nb::ndim<1>>>())
      .def("__getitem__", &test::test_xtadapt::operator[], "i"_a);

  nb::class_<emfields>(m, "emfields")
      .def("E", &emfields::E, "r"_a)
      .def("B", &emfields::B, "r"_a)
      .def("__repr__", &emfields::repr);

  nb::class_<emfields_uniform, emfields>(m, "emfields_uniform")
      .def(nb::init<double3, double3>(), "E_0"_a, "B_0"_a);

  nb::class_<emfields_dipole, emfields>(m, "emfields_dipole")
      .def(nb::init<double3>(), "m"_a);
}
