#include <gtest/gtest.h>
#include <tclap/CmdLine.h>
#include "../src/Random.hpp"
#include "../src/Data.hpp"

using namespace std;

TEST(ArithmeticTest, CorrectAddition) {
	EXPECT_EQ(2 + 2, 4);
}


// For the following tests, we know beforehand the data from the test files given on Moodle
// We test that each data has been well read and collected

TEST(DataReading, PopulationSize) {
	Data data("../data/test_input.txt","../data/test.fa");

	data.collectAll();

	EXPECT_EQ(data.getPopSize(), 21);
}

TEST(DataReading, NumberGenerations) {
	Data data("../data/test_input.txt","../data/test.fa");

	data.collectAll();

	EXPECT_EQ(data.getGenerations(), 3000);
}

TEST(DataReading, NumberAlleles) {
	Data data("../data/test_input.txt", "../data/test.fa");

	data.collectAll();

	EXPECT_EQ((int) data.getAllelesCount().size(), 2);
}

TEST(DataReading, MarkerSites) {
	Data data("../data/test_input.txt", "../data/test.fa");

	vector<int> knownMS = {1, 2, 3, 4};

	data.collectAll();

	EXPECT_EQ(data.getMarkerSites(), knownMS);
}


TEST(DataReading, InitialFrequencies) {
	Data data("../data/test_input.txt","../data/test.fa");

	vector<int> knownCount = {9, 12};

	data.collectAll();

	for (size_t i(0); i < (size_t) data.getAllelesCount().size(); ++i) {
		EXPECT_NEAR(data.getAllelesCount()[i], knownCount[i], 1E-3);
	}
}

TEST(DataReading, NumberReplicates) {
	Data data("../data/test_input.txt","../data/test.fa");

	data.collectAll();

	EXPECT_EQ(data.getReplicates(), 500);
}

TEST(DataReading, NucleotidesMutations) {
	Data data("../data/test_input.txt","../data/test.fa");

	vector <double> knownMut = {10E-8, 10E-8, 10E-8, 10E-8};

	data.collectAll();

	for (size_t i(0); i < data.getMutations().size() ; ++i) {
		EXPECT_NEAR(data.getMutations()[i], knownMut[i], 1E-3);
	}
}

TEST (DataReading, Bottleneck) {
	Data data("../data/test_input.txt","../data/test.fa");
	
	data.collectAll();
	
	EXPECT_EQ(data.getPopReduction(), 2);
	EXPECT_EQ(data.getBottleneckStart(), 20);
	EXPECT_EQ(data.getBottleneckEnd(), 40);
}


TEST(RandomTest, UniformDistribution) {
	double mean_uniform(0), input_mean(1.35), input_sd(2.8);

	RandomDist rng_unif(input_mean, input_sd, 10000, false);
	for (auto I : rng_unif.generate_numbers()) {
		EXPECT_GE(I, -3.5);
		EXPECT_LT(I, 6.2);
		mean_uniform += I*1e-4;
	}

	EXPECT_NEAR(input_mean, mean_uniform, 3 * input_sd / sqrt(1e4));
}

TEST(RandomTest, NormalDistribution) {
	double mean_normal(0), input_mean(1.35), input_sd(2.8);

	RandomDist rng_norm(input_mean, input_sd, 10000, true);
	for (auto I : rng_norm.generate_numbers()) {
		mean_normal += I * 1e-4;
	}

	EXPECT_NEAR(input_mean, mean_normal, 2 * input_sd / sqrt(1e4));
}





RandomDist* parse_args(int argc, char **argv) {
    TCLAP::CmdLine cmd("Random number generator");
    TCLAP::ValueArg<int> nsample("N", "sample_size", "Number of random numbers to generate",  true, 1, "int");
    cmd.add(nsample);
    TCLAP::ValueArg<double> meanval("m", "mean", "Mean value",  true, 0, "double");
    cmd.add( meanval );
    TCLAP::ValueArg< double > sdev("s", "std_dev", "Standard deviation",  true, 1, "double");
    cmd.add( sdev );
    TCLAP::SwitchArg  uniform("u", "uniform", "Uniform distribution", true);
    cmd.add( uniform );
    TCLAP::SwitchArg  normal("n", "normal", "Normal distribution", true);
    cmd.add(normal);
    cmd.parse( argc, argv );
    double m = meanval.getValue(), s = sdev.getValue();
    int ns = nsample.getValue();
    bool norm =  normal.getValue();
    if ( !(norm ^ uniform.getValue()) ) {
        throw(3);
    }
    RandomDist* rng = new RandomDist(m, s, ns, norm);
    return rng;
}

int main_random(int argc, char **argv) {
    RandomDist *rng = 0;
    try {
        rng = parse_args(argc, argv);
    } catch(const int n) {
        switch (n) {
        case 1:
            std::cerr << "Standard deviation must be positive\n";
            break;
        case 2:
            std::cerr << "Sample size must be positive\n";
            break;
        case 3:
            std::cerr << "Please choose either normal or uniform distribution\n";
            break;
        }
        return n;
    }
    for (auto I : rng->generate_numbers()) {
        std::cout << I << " "; std::cout << "\n";
    }
    if (rng) delete rng;
    return 0;
}


int main(int argc, char**argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
