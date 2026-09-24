
#pragma once

#include <openggcm/constants.hxx>
#include <openggcm/emfields.hxx>
#include <openggcm/tracing/particle.hxx>

#include <functional>
#include <string>

namespace openggcm
{

namespace tracing
{

class boris
{
public:
  boris(const emfields &emfields, double q, double m)
      : emfields_(emfields), q_(q), m_(m)
  {
  }

  std::string repr() const
  {
    return "boris(emfields=" + emfields_.get().repr() + ")";
  }

  void push_x(double3 &x, double3 &u, double dt) const
  {
    double inv_gamma = 1. / std::sqrt(1.0 + norm2(u));
    x += dt * u * inv_gamma * constants::c;
  }

  void push_u(double3 &x, double3 &u, double dt) const
  {
    auto E = emfields_.get().E(x);
    auto B = emfields_.get().B(x);

    double dq = 0.5 * (q_ / m_) * dt;
    // Half acceleration due to electric field
    double3 um = u + dq * E / constants::c;

    // Rotation due to magnetic field
    double root = dq / sqrt(1. + sqr(um[0]) + sqr(um[1]) + sqr(um[2]));
    double3 tau = root * B;
    double tau_norm = 1. / (1. + sqr(tau[0]) + sqr(tau[1]) + sqr(tau[2]));
    double3 up;
    up[0] = ((1. + sqr(tau[0]) - sqr(tau[1]) - sqr(tau[2])) * um[0] +
             (2. * tau[0] * tau[1] + 2. * tau[2]) * um[1] +
             (2. * tau[0] * tau[2] - 2. * tau[1]) * um[2]) *
            tau_norm;
    up[1] = ((2. * tau[0] * tau[1] - 2. * tau[2]) * um[0] +
             (1. - sqr(tau[0]) + sqr(tau[1]) - sqr(tau[2])) * um[1] +
             (2. * tau[1] * tau[2] + 2. * tau[0]) * um[2]) *
            tau_norm;
    up[2] = ((2. * tau[0] * tau[2] + 2. * tau[1]) * um[0] +
             (2. * tau[1] * tau[2] - 2. * tau[0]) * um[1] +
             (1. - sqr(tau[0]) - sqr(tau[1]) + sqr(tau[2])) * um[2]) *
            tau_norm;

    u = up + dq * E / constants::c;
  }

  void push(particles &prts, std::optional<double> t_final,
            std::optional<int> max_steps, std::optional<double> dt_max,
            double dt_max_gyro) const
  {
    double qprime = 0.5 * q_ / m_;

    if (!t_final.has_value() && !max_steps.has_value())
    {
      throw std::runtime_error(
          "boris::push: either t_final or max_steps must be specified");
    }

    for (int n = 0; n < prts.size(); ++n)
    {
      int step = 0;
      while (true)
      {
        double &t = prts.t(n);
        double3 &x = prts.r(n);
        double3 &u = prts.u(n);

        double3 B = emfields_.get().B(x);
        double gamma = std::sqrt(1.0 + norm2(u));
        double om_c = 2.0 * std::abs(qprime) * norm(B) / gamma;
        double dt = dt_max_gyro * 2.0 * constants::pi / om_c;
        if (dt_max.has_value())
        {
          dt = std::min(dt_max.value(), dt);
        }

        if (t_final.has_value() && t >= t_final.value())
        {
          break;
        }
        if (max_steps.has_value() && step >= max_steps.value())
        {
          break;
        }

        push_x(x, u, .5 * dt);
        push_u(x, u, dt);
        push_x(x, u, .5 * dt);
        t += dt;
        step++;
      }
    }
  }

private:
  std::reference_wrapper<const emfields>
      emfields_; // FIXME potential lifetime issues with nanobind
  double q_, m_;
};

} // namespace tracing
} // namespace openggcm
