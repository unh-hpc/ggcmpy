
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
  particle(double3 x, double3 u) : x(x), u(u) {}

  double3 v() const { return (constants::c / gamma()) * u; }

  double gamma() const
  {
    return std::sqrt(1.0 + sqr(u[0]) + sqr(u[1]) + sqr(u[2]));
  }

  std::string repr() const
  {
    return "particle(x=[" + std::to_string(x[0]) + ", " + std::to_string(x[1]) +
           ", " + std::to_string(x[2]) + "], u=[" + std::to_string(u[0]) +
           ", " + std::to_string(u[1]) + ", " + std::to_string(u[2]) + "])";
  }

  double3 x;
  double3 u;
};

class particles
{
public:
  void add_particle(double t, double3 r, double3 u)
  {
    t_.push_back(t);
    r_.push_back(r);
    u_.push_back(u);
  }

  std::tuple<std::vector<double>, std::vector<double3>, std::vector<double3>>
  to_tuple() const
  {
    return std::make_tuple(t_, r_, u_);
  }

private:
  std::vector<double> t_;
  std::vector<double3> r_;
  std::vector<double3> u_;
};

} // namespace tracing
} // namespace openggcm
