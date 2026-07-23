
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

template <typename Container> class cic_weights
{
public:
  cic_weights(const Container &crd) : crd_(crd) {}

  std::pair<std::size_t, double2> operator()(double x) const
  {
    auto it = std::upper_bound(crd_.data(), crd_.data() + crd_.size(), x);
    std::size_t i = std::distance(crd_.data(), it) - 1;
    if (i >= crd_.size() - 1)
    {
      return {std::size_t(-1), {NaN, NaN}};
    }
    double w1 = (x - crd_(i)) / (crd_(i + 1) - crd_(i));
    double w0 = 1.0 - w1;
    return {i, {w0, w1}};
  }

private:
  const Container &crd_;
};

} // namespace interpolate

template <typename Container1d, typename Container3d>
class emfields_yee : public emfields
{
public:
  using array1d = Container1d;
  using array3d = Container3d;

  emfields_yee(array3d e1x, array3d e1y, array3d e1z, array3d b1x, array3d b1y,
               array3d b1z, array1d x, array1d y, array1d z, array1d x_nc,
               array1d y_nc, array1d z_nc)
      : e1x_(e1x), e1y_(e1y), e1z_(e1z), b1x_(b1x), b1y_(b1y), b1z_(b1z), x_(x),
        y_(y), z_(z), x_nc_(x_nc), y_nc_(y_nc), z_nc_(z_nc)
  {
  }

protected:
  array3d e1x_, e1y_, e1z_;
  array3d b1x_, b1y_, b1z_;
  array1d x_, y_, z_;
  array1d x_nc_, y_nc_, z_nc_;
};

template <typename Container1d, typename Container3d>
class emfields_yee_cic : public emfields_yee<Container1d, Container3d>
{
public:
  using base_type = emfields_yee<Container1d, Container3d>;
  using typename base_type::array1d;
  using typename base_type::array3d;

  using base_type::base_type;

  double3 E(double3 r) const override
  {
    return {interpolate(r, this->e1x_, this->x_, this->y_nc_, this->z_nc_),
            interpolate(r, this->e1y_, this->x_nc_, this->y_, this->z_nc_),
            interpolate(r, this->e1z_, this->x_nc_, this->y_nc_, this->z_)};
  }

  double3 B(double3 r) const override
  {
    return {interpolate(r, this->b1x_, this->x_nc_, this->y_, this->z_),
            interpolate(r, this->b1y_, this->x_, this->y_nc_, this->z_),
            interpolate(r, this->b1z_, this->x_, this->y_, this->z_nc_)};
  }

  std::string repr() const override { return "emfields_yee_cic()"; }

private:
  static double interpolate(const double3 &r, const array3d &f,
                            const array1d &x, const array1d &y,
                            const array1d &z)
  {
    auto [ix, wx] = interpolate::cic_weights(x)(r(0));
    auto [iy, wy] = interpolate::cic_weights(y)(r(1));
    auto [iz, wz] = interpolate::cic_weights(z)(r(2));

    return (f(ix + 0, iy + 0, iz + 0) * wx[0] * wy[0] * wz[0] +
            f(ix + 1, iy + 0, iz + 0) * wx[1] * wy[0] * wz[0] +
            f(ix + 0, iy + 1, iz + 0) * wx[0] * wy[1] * wz[0] +
            f(ix + 1, iy + 1, iz + 0) * wx[1] * wy[1] * wz[0] +
            f(ix + 0, iy + 0, iz + 1) * wx[0] * wy[0] * wz[1] +
            f(ix + 1, iy + 0, iz + 1) * wx[1] * wy[0] * wz[1] +
            f(ix + 0, iy + 1, iz + 1) * wx[0] * wy[1] * wz[1] +
            f(ix + 1, iy + 1, iz + 1) * wx[1] * wy[1] * wz[1]);
  }
};

} // namespace openggcm
