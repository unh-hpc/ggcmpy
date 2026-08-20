
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
  particles() = default;
  template <typename E>
  particles(const xt::xexpression<E> &t, const xt::xexpression<E> &r,
            const xt::xexpression<E> &u)
  {
    t_.reserve(t.derived_cast().shape(0));
    r_.reserve(r.derived_cast().shape(0));
    u_.reserve(u.derived_cast().shape(0));
    for (std::size_t i = 0; i < t.derived_cast().shape(0); ++i)
    {
      t_.push_back(t.derived_cast()(i));
      r_.push_back(double3{r.derived_cast()(i, 0), r.derived_cast()(i, 1),
                           r.derived_cast()(i, 2)});
      u_.push_back(double3{u.derived_cast()(i, 0), u.derived_cast()(i, 1),
                           u.derived_cast()(i, 2)});
    }
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
