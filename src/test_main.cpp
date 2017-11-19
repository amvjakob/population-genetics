#include "gtest/gtest.h"
#include "random.hpp"
#include <cmath>

TEST(ArithmeticTest, CorrectAddition) {
	EXPECT_EQ(2 + 2, 4);
}

TEST(randomTest, randomDistributions) {

double mean_uniform(0), mean_normal(0), input_mean(1.35), input_sd(2.8);
Randomdist rng_unif(input_mean, input_sd, 10000, false, 12345);

for (auto I : rng_unif.generate_numbers()) {
EXPECT_GE(I, -3.5);
EXPECT_LT(I, 6.2);
mean_uniform += I*1e-4;
}

EXPECT_NEAR(input_mean, mean_uniform, 3*input_sd / sqrt(1e4));
Randomdist rng_norm(input_mean, input_sd, 10000, true, 12345);
for (auto I : rng_norm.generate_numbers()) mean_normal += I*1e-4;
EXPECT_NEAR(input_mean, mean_normal, 2*input_sd / sqrt(1e4));
}

int main(int argc, char **argv) {
::testing::InitGoogleTest(&argc, argv);
return RUN_ALL_TESTS();
}

