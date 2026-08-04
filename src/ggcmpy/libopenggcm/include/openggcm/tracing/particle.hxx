
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

  double t(std::size_t i) const { return t_[i]; }
  double3 r(std::size_t i) const { return r_[i]; }
  double3 u(std::size_t i) const { return u_[i]; }
  double &t(std::size_t i) { return t_[i]; }
  double3 &r(std::size_t i) { return r_[i]; }
  double3 &u(std::size_t i) { return u_[i]; }

  std::vector<double> &t() { return t_; }
  std::vector<double3> &r() { return r_; }
  std::vector<double3> &u() { return u_; }

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
