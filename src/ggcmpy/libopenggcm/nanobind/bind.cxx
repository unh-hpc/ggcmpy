
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include <openggcm/emfields.hxx>

#include <nanobind/stl/detail/nb_array.h>
#include <xtensor/containers/xadapt.hpp>
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
  using nb_array_type = nb::ndarray<double, nb::c_contig>;
  using xt_adapter_type = decltype(xt::adapt_smart_ptr(
      std::declval<double *>(), std::declval<std::vector<std::size_t>>(),
      std::declval<std::unique_ptr<nb_array_type>>()));

  std::vector<std::size_t> make_shape(const nb_array_type &a)
  {
    std::vector<std::size_t> shape(a.ndim());
    for (std::size_t i = 0; i < a.ndim(); ++i)
    {
      shape[i] = a.shape(i);
    }
    return shape;
  }

public:
  test_xtadapt(nb_array_type a)
      : a_xt_(std::move(xt::adapt_smart_ptr(
            a.data(), make_shape(a), std::make_unique<nb_array_type>(a))))
  {
  }

  double operator()(std::size_t i) const { return a_xt_(i); }
  double operator()(std::size_t i, std::size_t j) const { return a_xt_(i, j); }

  std::string repr() const
  {
    std::string rv = "test_xtadapt(shape=[";
    auto shape = a_xt_.shape();
    for (std::size_t i = 0; i < shape.size(); ++i)
    {
      rv += std::to_string(shape[i]);
      if (i < shape.size() - 1)
      {
        rv += ", ";
      }
    }
    rv += "])";
    return rv;
  }

private:
  xt_adapter_type a_xt_;
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
      .def(nb::init<nb::ndarray<double, nb::c_contig>>(), nb::arg().noconvert())
      .def("__getitem__",
           [](test::test_xtadapt &self, std::size_t i) { return self(i); })
      .def("__getitem__",
           [](test::test_xtadapt &self, std::tuple<std::size_t> idx)
           { return self(std::get<0>(idx)); })
      .def("__getitem__", [](test::test_xtadapt &self,
                             std::tuple<std::size_t, std::size_t> idx)
           { return self(std::get<0>(idx), std::get<1>(idx)); })
      .def("__repr__", &test::test_xtadapt::repr);

  nb::class_<emfields>(m, "emfields")
      .def("E", &emfields::E, "r"_a)
      .def("B", &emfields::B, "r"_a)
      .def("__repr__", &emfields::repr);

  nb::class_<emfields_uniform, emfields>(m, "emfields_uniform")
      .def(nb::init<double3, double3>(), "E_0"_a, "B_0"_a);

  nb::class_<emfields_dipole, emfields>(m, "emfields_dipole")
      .def(nb::init<double3>(), "m"_a);
}
