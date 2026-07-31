
#pragma once

#include <openggcm/constants.hxx>
#include <openggcm/util.hxx>

namespace openggcm
{
namespace tracing
{

class particle
{
public:
  particle(double3 x, double3 u, double q, double m) : x(x), u(u), q(q), m(m) {}

  double3 v() const { return (constants::c / gamma()) * u; }

  double gamma() const
  {
    return std::sqrt(1.0 + sqr(u[0]) + sqr(u[1]) + sqr(u[2]));
  }

  std::string repr() const
  {
    return "particle(x=[" + std::to_string(x[0]) + ", " + std::to_string(x[1]) +
           ", " + std::to_string(x[2]) + "], u=[" + std::to_string(u[0]) +
           ", " + std::to_string(u[1]) + ", " + std::to_string(u[2]) +
           "], q=" + std::to_string(q) + ", m=" + std::to_string(m) + ")";
  }

  double3 x;
  double3 u;
  double q;
  double m;
};

} // namespace tracing
} // namespace openggcm
