
#pragma once

#include <array>

namespace openggcm
{

class emfields_constant
{
public:
  emfields_constant(std::array<double, 3> E0, std::array<double, 3> B0)
    : _E0(E0), _B0(B0)
  {
  }

  std::array<double, 3> E(double x, double y, double z) const { return _E0; }
  std::array<double, 3> B(double x, double y, double z) const { return _B0; }

private:
  std::array<double, 3> _E0;
  std::array<double, 3> _B0;
};

} // namespace openggcm
