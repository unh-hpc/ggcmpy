
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

  void push_x(particle &p, double dt) const { p.x += dt * p.v(); }

  void push_v(particle &p, double dt) const
  {
    auto E = emfields_.get().E(p.x);
    auto B = emfields_.get().B(p.x);

    double dq = 0.5 * (q_ / m_) * dt;
    // Half acceleration due to electric field
    double3 um = p.u + dq * E / constants::c;

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

    p.u = up + dq * E / constants::c;
  }

  void push(particle &prt, double t_max, double dt_max, double gyro_max,
            particles &prts) const
  {
    double qprime = 0.5 * q_ / m_;
    double3 B = emfields_.get().B(prt.x);
    double om_c = 2.0 * std::abs(qprime) * norm(B);
    double dt = std::min(dt_max, gyro_max * 2.0 * constants::pi / om_c);

    double t = 0.;
    prts.add_particle(t, prt.x, prt.u);

    while (t < t_max)
    {
      push_x(prt, .5 * dt);
      push_v(prt, dt);
      push_x(prt, .5 * dt);
      t += dt;
      prts.add_particle(t, prt.x, prt.u);
    }
  }

private:
  std::reference_wrapper<const emfields>
      emfields_; // FIXME potential lifetime issues with nanobind
  double q_, m_;
};

} // namespace tracing
} // namespace openggcm
