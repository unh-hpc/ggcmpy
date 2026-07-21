
#pragma once

#include <numbers>

namespace openggcm
{
namespace constants
{

constexpr double pi = std::numbers::pi;
constexpr double c = 299792458.0;       // speed of light in m/s
constexpr double e = 1.602176634e-19;   // elementary charge in C
constexpr double m_e = 9.10938356e-31;  // electron mass in kg
constexpr double mu0 = 4.0 * pi * 1e-7; // vacuum permeability

} // namespace constants
} // namespace openggcm
