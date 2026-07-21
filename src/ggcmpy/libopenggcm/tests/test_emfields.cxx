
#include <gtest/gtest.h>

#include <openggcm/emfields.hxx>

TEST(EmfieldsTest, ConstantFields)
{
  double3 E_0 = {1.0, 2.0, 3.0};
  double3 B_0 = {4.0, 5.0, 6.0};
  openggcm::emfields_uniform fields(E_0, B_0);

  EXPECT_EQ(fields.E(0.0, 0.0, 0.0), E_0);
  EXPECT_EQ(fields.B(0.0, 0.0, 0.0), B_0);
}
