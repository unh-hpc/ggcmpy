
#pragma once

#include <openggcm/constants.hxx>
#include <openggcm/util.hxx>

#include <xtensor/containers/xtensor.hpp>

#include <string>

namespace openggcm
{

class emfields
{
public:
  virtual ~emfields() = default;

  virtual double3 E(double3 r) const = 0;
  virtual double3 B(double3 r) const = 0;
  virtual std::string repr() const = 0;
};

class emfields_uniform : public emfields
{
public:
  emfields_uniform(double3 E_0, double3 B_0) : _E_0(E_0), _B_0(B_0) {}

  double3 E(double3 r) const override { return _E_0; }
  double3 B(double3 r) const override { return _B_0; }

  std::string repr() const override
  {
    return "emfields_uniform(E_0=[" + std::to_string(_E_0[0]) + ", " +
           std::to_string(_E_0[1]) + ", " + std::to_string(_E_0[2]) +
           "], B_0=[" + std::to_string(_B_0[0]) + ", " +
           std::to_string(_B_0[1]) + ", " + std::to_string(_B_0[2]) + "])";
  }

private:
  double3 _E_0;
  double3 _B_0;
};

class emfields_dipole : public emfields
{
public:
  emfields_dipole(double3 m) : m_(m) {}

  double3 E(double3 r) const override { return {}; }
  double3 B(double3 r) const override
  {
    double rnorm = norm(r);
    double3 rhat = r / rnorm;
    return constants::mu0 / (4.0 * constants::pi) *
           (3.0 * dot(m_, rhat) * rhat - m_) / std::pow(rnorm, 3);
  }

  std::string repr() const override
  {
    return "emfields_dipole(m=[" + std::to_string(m_[0]) + ", " +
           std::to_string(m_[1]) + ", " + std::to_string(m_[2]) + "])";
  }

private:
  double3 m_; // magnetic moment
};

namespace interpolate
{

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

class cic_weights
{
public:
  cic_weights(const xt::xtensor<double, 1> &crd) : crd_(crd) {}

  std::pair<std::size_t, double2> operator()(double x) const
  {
    auto it = std::upper_bound(crd_.cbegin(), crd_.cend(), x);
    std::size_t i = std::distance(crd_.cbegin(), it) - 1;
    if (i >= crd_.size() - 1)
    {
      return {std::size_t(-1), {NaN, NaN}};
    }
    double w1 = (x - crd_(i)) / (crd_(i + 1) - crd_(i));
    double w0 = 1.0 - w1;
    return {i, {w0, w1}};
  }

private:
  xt::xtensor<double, 1> crd_;
};

} // namespace interpolate

} // namespace openggcm
