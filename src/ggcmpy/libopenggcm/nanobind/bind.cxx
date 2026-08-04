
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <openggcm/emfields.hxx>
#include <openggcm/tracing/boris.hxx>
#include <openggcm/tracing/particle.hxx>

#include <nanobind/stl/detail/nb_array.h>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xfixed.hpp>

template <typename T>
using nb_ndarray_type = nanobind::ndarray<T, nanobind::any_contig>;

namespace detail
{
template <typename T>
std::vector<std::size_t> shape_from_nb(const nb_ndarray_type<T> &a)
{
  std::vector<std::size_t> shape(a.ndim());
  for (std::size_t i = 0; i < a.ndim(); ++i)
  {
    shape[i] = a.shape(i);
  }
  return shape;
}

template <typename T>
xt::layout_type layout_from_nb(const nb_ndarray_type<T> &a)
{
  size_t ndim = a.ndim();

  // 0D or 1D arrays are contiguous in both definitions
  if (ndim <= 1)
  {
    return xt::layout_type::row_major; // F as well...
  }

  // Expected stride for the next dimension
  std::size_t expected_c_stride = 1;
  std::size_t expected_f_stride = 1;

  // Check C-Contiguous (Row-Major): Strides increase from right to left
  bool is_c_contig = true;
  for (int i = (int)ndim - 1; i >= 0; --i)
  {
    if (a.stride(i) != expected_c_stride)
    {
      is_c_contig = false;
    }
    expected_c_stride *= a.shape(i);
  }
  if (is_c_contig)
  {
    return xt::layout_type::row_major;
  }

  // Check F-Contiguous (Column-Major): Strides increase from left to right
  bool is_f_contig = true;
  for (size_t i = 0; i < ndim; ++i)
  {
    if (a.stride(i) != expected_f_stride)
    {
      is_f_contig = false;
    }
    expected_f_stride *= a.shape(i);
  }
  if (is_f_contig)
  {
    return xt::layout_type::column_major;
  }

  return xt::layout_type::any; // should never happen
}

} // namespace detail

template <typename T> auto xt_adapt_ndarray(nb_ndarray_type<T> arr)
{
  return xt::adapt_smart_ptr<xt::layout_type::dynamic>(
      arr.data(), detail::shape_from_nb(arr),
      std::make_unique<nb_ndarray_type<T>>(arr), detail::layout_from_nb(arr));
}

// template <typename T>
// using xt_ndarray =
//   decltype(xt_adapt_ndarray(std::declval<nb_ndarray_type<T>>()));

template <typename T>
using xt_ndarray =
    xt::xarray_adaptor<xt::xbuffer_adaptor<T *, xt::smart_ownership,
                                           std::unique_ptr<nb_ndarray_type<T>>>,
                       xt::layout_type::dynamic, std::vector<std::size_t>>;

NAMESPACE_BEGIN(NB_NAMESPACE)
NAMESPACE_BEGIN(detail)

template <typename Type, size_t Size>
struct type_caster<xt::xtensor_fixed<Type, xt::xshape<Size>>>
    : array_caster<xt::xtensor_fixed<Type, xt::xshape<Size>>, Type, Size>
{
};

template <typename T> struct type_caster<xt::xarray<T>>
{
  NB_TYPE_CASTER(xt::xarray<T>,
                 const_name("xt::xarray<") + const_name<T>() + const_name(">"));

  bool from_python(handle src, uint8_t flags, cleanup_list *cleanup) noexcept
  {
    auto ndarray_caster = make_caster<nb_ndarray_type<T>>();
    if (!ndarray_caster.from_python(src, flags, cleanup))
    {
      return false;
    }

    value = xt_adapt_ndarray(ndarray_caster.value);
    return true;
  };
};

// template <typename T>
// struct type_caster<xt_ndarray<T>>
// {
//   NB_TYPE_CASTER(xt_ndarray<T>, const_name("xt_ndarray<") + const_name<T>() +
//                                   const_name(">&&"));

//   bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept
//   {
//     auto ndarray_caster = make_caster<nb_ndarray_type<T>>();
//     if (!ndarray_caster.from_python(src, flags, cleanup)) {
//       return false;
//     }

//     value = xt_adapt_ndarray(ndarray_caster.value);
//     return true;
//   };
// };

NAMESPACE_END(detail)
NAMESPACE_END(NB_NAMESPACE)

namespace nb = nanobind;
using namespace nb::literals;
using namespace openggcm;
using namespace openggcm::tracing;

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
  test_xtadapt(nb_ndarray_type<double> a) : a_xt_(xt_adapt_ndarray(a)) {}

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
  xt_ndarray<double> a_xt_;
};
} // namespace test

NB_MODULE(_openggcm, m)
{
  m.def("xt_fixed_from_python",
        [](xt::xtensor_fixed<double, xt::xshape<3>> a) { return a; });

  m.def("xt_array_from_python",
        [](xt::xarray<double> a)
        {
          return "xt::array ndim=" + std::to_string(a.dimension());
          a(1) = 99.;
        });

  m.def("xt_ndarray_from_python_manual",
        [](nb_ndarray_type<double> a)
        {
          auto a_xt = xt_adapt_ndarray(a);
          a_xt(1) = 99.;
          return "xt_ndarray ndim=" + std::to_string(a_xt.dimension());
        });

  // m.def("xt_ndarray_from_python", [](xt_ndarray<double>&& a) {
  //   a(1) = 99.;
  //   return "xt_ndarray ndim=" + std::to_string(a.dimension());
  // });

  nb::class_<test::test_ndarray>(m, "test_ndarray")
      .def(nb::init<nb::ndarray<double, nb::ndim<1>>>())
      .def("__getitem__", &test::test_ndarray::operator[], "i"_a);

  nb::class_<test::test_xtadapt>(m, "test_xtadapt")
      .def(nb::init<nb_ndarray_type<double>>(), nb::arg().noconvert())
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
      .def(nb::init<double3, double3>(), "E_0"_a = double3{0.0, 0.0, 0.0},
           "B_0"_a = double3{0.0, 0.0, 0.0});

  nb::class_<emfields_dipole, emfields>(m, "emfields_dipole")
      .def(nb::init<double3>(), "m"_a);

  using array1d = nb::ndarray<const double, nb::ndim<1>, nb::any_contig>;
  using array3d = nb::ndarray<const double, nb::ndim<3>, nb::any_contig>;

  nb::class_<emfields_yee<array1d, array3d>, emfields>(m, "emfields_yee");

  nb::class_<emfields_yee_cic<array1d, array3d>,
             emfields_yee<array1d, array3d>>(m, "emfields_yee_cic")
      .def(nb::init<array3d, array3d, array3d, array3d, array3d, array3d,
                    array1d, array1d, array1d, array1d, array1d, array1d>(),
           "e1x"_a, "e1y"_a, "e1z"_a, "b1x"_a, "b1y"_a, "b1z"_a, "x"_a, "y"_a,
           "z"_a, "x_nc"_a, "y_nc"_a, "z_nc"_a);

  nb::class_<emfields_yee_tsc<array1d, array3d>,
             emfields_yee<array1d, array3d>>(m, "emfields_yee_tsc")
      .def(nb::init<array3d, array3d, array3d, array3d, array3d, array3d,
                    array1d, array1d, array1d, array1d, array1d, array1d>(),
           "e1x"_a, "e1y"_a, "e1z"_a, "b1x"_a, "b1y"_a, "b1z"_a, "x"_a, "y"_a,
           "z"_a, "x_nc"_a, "y_nc"_a, "z_nc"_a);

  // ==================================================================
  // tracing submodule
  nb::module_ tracing = m.def_submodule("tracing", "tracing submodule");

  // ------------------------------------------------------------------
  // particle
  nb::class_<particle>(tracing, "particle")
      .def(nb::init<double3, double3, double, double>(), "x"_a, "u"_a, "q"_a,
           "m"_a)
      .def_prop_ro("v", [](const particle &p) { return p.v(); })
      .def_prop_ro("gamma", [](const particle &p) { return p.gamma(); })
      .def("__repr__", &particle::repr)
      .def_rw("x", &particle::x)
      .def_rw("u", &particle::u)
      .def_rw("q", &particle::q)
      .def_rw("m", &particle::m);

  // ------------------------------------------------------------------
  // particles
  nb::class_<particles>(tracing, "particles")
      .def(nb::init<>())
      .def("add_particle", &particles::add_particle, "t"_a, "r"_a, "u"_a)
      .def("to_tuple", &particles::to_tuple);

  // ------------------------------------------------------------------
  // boris
  nb::class_<boris>(tracing, "boris")
      .def(nb::init<const emfields &>(), "emfields"_a)
      .def("__repr__", &boris::repr)
      .def("push", &boris::push, "t"_a, "p"_a, "prts"_a, "dt_max"_a,
           "gyro_max"_a);
}
