#include "../src/data.hpp"
#include "googletest/include/gtest/gtest.h"

using namespace std;

//For the following tests, we know beforehand the data from the test files given on Moodle
// We test that each data has been well read and collected

TEST(DataReading, PopulationSize) {
	Data data("../data/input.txt","../data/test.fa");
	
	data.collectAll();
	
	EXPECT_EQ(data.getPopSize(), 21);
}

TEST(DataReading, NumberGenerations) {
	
	Data data("../data/input.txt","../data/test.fa");
	
	data.collectAll();
	
	EXPECT_EQ(data.getGenerations(), 3000);
	
}

TEST(DataReading, NumberAlleles) {
	
	Data data("../data/input.txt","../data/test.fa");
	
	data.collectAll();
	
	EXPECT_EQ(data.getNumberAlleles(), 2);
	
}

TEST(DataReading, MarkerSites) {
	
	Data data("../data/input.txt","../data/test.fa");
	
	vector <double> knownMS = {1,2,3,4};
	
	data.collectAll();
	
	EXPECT_EQ(data.getMarkerSites(), knownMS);
	
}


TEST(DataReading, InitialFrequencies) {
	
	Data data("../data/input.txt","../data/test.fa");
	
	vector <double> knownFq = {0.428571, 0.571429};
	
	data.collectAll();
	
	for (size_t i(0); i < data.getNumberAlleles(); ++i) {
	EXPECT_NEAR(data.getAllelesFq()[i], knownFq[i], 1E-3);
	}
	
}

TEST(DataReading, NumberReplicates) {
	
	Data data("../data/input.txt","../data/test.fa");
	
	data.collectAll();
	
	EXPECT_EQ(data.getReplicates(), 500);
	
}

TEST(DataReading, NucleotidesMutations) {
	
	Data data("../data/input.txt","../data/test.fa");
	
	vector <double> knownMut = {10E-8, 10E-8, 10E-8, 10E-8};
	
	data.collectAll();
	
	for (size_t i(0); i < data.getMutations().size() ; ++i) {
	EXPECT_NEAR(data.getMutations()[i], knownMut[i], 1E-3);
	}
	
}

int main(int argc, char**argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
