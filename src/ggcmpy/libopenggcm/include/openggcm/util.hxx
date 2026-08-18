
#pragma once

#include <xtensor/containers/xfixed.hpp>

namespace openggcm
{

using double2 = xt::xtensor_fixed<double, xt::xshape<2>>;
using double3 = xt::xtensor_fixed<double, xt::xshape<3>>;

inline double sqr(double x) { return x * x; }

inline double norm(double3 x)
{
  return std::sqrt(sqr(x[0]) + sqr(x[1]) + sqr(x[2]));
}

inline double dot(double3 x, double3 y)
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

} // namespace openggcm
