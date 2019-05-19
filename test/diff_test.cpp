#include "differentiation.hpp"

#include "gtest/gtest.h"

TEST(TestDiff1D, Pattern_f_x)
{
    double f_x = numcalc::diff<numcalc::f_x>(exp, 0.0);
    EXPECT_NEAR(f_x, 1.0, 1e-5);
}

TEST(TestDiff1D, Pattern_f_xx)
{
    double f_xx = numcalc::diff<numcalc::f_xx>(log, 2.0);
    EXPECT_NEAR(f_xx, -0.25, 1e-5);
}

TEST(TestDiff2D, Pattern_f_x)
{
    auto f = [](double x, double y){return exp(x + 2.0*y);};
    double f_x = numcalc::diff<numcalc::f_x>(f, 1.0, 1.0);
    EXPECT_NEAR(f_x, exp(3.0), 1e-5);
}

TEST(TestDiff2D, Pattern_f_y)
{
    auto f = [](double x, double y){return exp(x + 2.0*y);};
    double f_y = numcalc::diff<numcalc::f_y>(f, 1.0, 1.0);
    EXPECT_NEAR(f_y, 2.0*exp(3.0), 1e-5);
}

TEST(TestDiff2D, Pattern_f_xx)
{
    auto f = [](double x, double y){return exp(x + 2.0*y);};
    double f_xx = numcalc::diff<numcalc::f_xx>(f, 1.0, 1.0);
    EXPECT_NEAR(f_xx, exp(3.0), 1e-5);
}

TEST(TestDiff2D, Pattern_f_xy)
{
    auto f = [](double x, double y){return exp(x + 2.0*y);};
    double f_xy = numcalc::diff<numcalc::f_xy>(f, 1.0, 1.0);
    EXPECT_NEAR(f_xy, 2.0*exp(3.0), 1e-5);
}

TEST(TestDiff2D, Pattern_f_yy)
{
    auto f = [](double x, double y){return exp(x + 2.0*y);};
    double f_yy = numcalc::diff<numcalc::f_yy>(f, 1.0, 1.0);
    EXPECT_NEAR(f_yy, 4.0*exp(3.0), 1e-5);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
