
#pragma once

#include <array>

namespace openggcm
{

class emfields_uniform
{
public:
  emfields_uniform(std::array<double, 3> E_0, std::array<double, 3> B_0)
    : _E_0(E_0), _B_0(B_0)
  {
  }

  std::array<double, 3> E(double x, double y, double z) const { return _E_0; }
  std::array<double, 3> B(double x, double y, double z) const { return _B_0; }

private:
  std::array<double, 3> _E_0;
  std::array<double, 3> _B_0;
};

} // namespace openggcm
