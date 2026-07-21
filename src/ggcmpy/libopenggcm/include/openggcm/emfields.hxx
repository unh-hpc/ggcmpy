
#pragma once

#include <openggcm/util.hxx>

#include <string>

namespace openggcm
{

class emfields_uniform
{
public:
  emfields_uniform(double3 E_0, double3 B_0) : _E_0(E_0), _B_0(B_0) {}

  double3 E(double x, double y, double z) const { return _E_0; }
  double3 B(double x, double y, double z) const { return _B_0; }

  std::string repr() const
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

} // namespace openggcm
