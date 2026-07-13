
#include <gtest/gtest.h>

#include <openggcm/emfields.hxx>

TEST(EmfieldsTest, ConstantFields)
{
  std::array<double, 3> E0 = {1.0, 2.0, 3.0};
  std::array<double, 3> B0 = {4.0, 5.0, 6.0};
  openggcm::emfields_constant fields(E0, B0);

  EXPECT_EQ(fields.E(0.0, 0.0, 0.0), E0);
  EXPECT_EQ(fields.B(0.0, 0.0, 0.0), B0);
}
