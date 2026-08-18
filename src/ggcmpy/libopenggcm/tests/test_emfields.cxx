
#include <gtest/gtest.h>

#include <openggcm/emfields.hxx>

using namespace openggcm;

TEST(EmfieldsTest, ConstantFields)
{
  double3 E_0 = {1.0, 2.0, 3.0};
  double3 B_0 = {4.0, 5.0, 6.0};
  emfields_uniform fields(E_0, B_0);

  EXPECT_EQ(fields.E({0.0, 0.0, 0.0}), E_0);
  EXPECT_EQ(fields.B({0.0, 0.0, 0.0}), B_0);
}

TEST(interpolate, cic_weights)
{
  xt::xtensor<double, 1> crd = {1.0, 2.0, 4.0};
  interpolate::cic_weights weights(crd);

  {
    auto [i, w] = weights(2.5);
    EXPECT_EQ(i, 1);
    EXPECT_NEAR(w[0], 0.75, 1e-10);
    EXPECT_NEAR(w[1], 0.25, 1e-10);
  }

  {
    auto [i, w] = weights(.5);
    EXPECT_EQ(i, std::size_t(-1));
    EXPECT_TRUE(std::isnan(w[0]));
    EXPECT_TRUE(std::isnan(w[1]));
  }

  {
    auto [i, w] = weights(4.5);
    EXPECT_EQ(i, std::size_t(-1));
    EXPECT_TRUE(std::isnan(w[0]));
    EXPECT_TRUE(std::isnan(w[1]));
  }
}

TEST(interpolate, tsc_weights)
{
  xt::xtensor<double, 1> crd = {1.0, 2.0, 3.0};
  interpolate::tsc_weights weights(crd);

  {
    auto [i, w] = weights(2.0);
    EXPECT_EQ(i, 0);
    EXPECT_NEAR(w[0], 0.125, 1e-10);
    EXPECT_NEAR(w[1], 0.75, 1e-10);
    EXPECT_NEAR(w[2], 0.125, 1e-10);
  }

  {
    auto [i, w] = weights(1.5);
    EXPECT_EQ(i, 0);
    EXPECT_NEAR(w[0], 0.5, 1e-10);
    EXPECT_NEAR(w[1], 0.5, 1e-10);
    EXPECT_NEAR(w[2], 0.0, 1e-10);
  }

  {
    auto [i, w] = weights(1.2);
    EXPECT_EQ(i, std::size_t(-1));
    EXPECT_TRUE(std::isnan(w[0]));
    EXPECT_TRUE(std::isnan(w[1]));
    EXPECT_TRUE(std::isnan(w[2]));
  }

  {
    auto [i, w] = weights(2.8);
    EXPECT_EQ(i, std::size_t(-1));
    EXPECT_TRUE(std::isnan(w[0]));
    EXPECT_TRUE(std::isnan(w[1]));
    EXPECT_TRUE(std::isnan(w[2]));
  }
}
