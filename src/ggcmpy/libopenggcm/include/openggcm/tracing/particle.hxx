
#pragma once

#include <openggcm/constants.hxx>
#include <openggcm/util.hxx>

namespace openggcm
{
namespace tracing
{

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
