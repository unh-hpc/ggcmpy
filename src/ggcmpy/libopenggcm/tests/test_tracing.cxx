
#include <gtest/gtest.h>

#include <openggcm/tracing/particle.hxx>

using namespace openggcm;

TEST(particle, ctor)
{
  tracing::particle p({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0});
  EXPECT_EQ(p.x, double3({1.0, 2.0, 3.0}));
  EXPECT_EQ(p.u, double3({4.0, 5.0, 6.0}));
}

TEST(particle, gamma)
{
  tracing::particle p({1.0, 2.0, 3.0}, {0.0, 0.0, 1.0});
  EXPECT_EQ(p.gamma(), std::sqrt(2.));
}
