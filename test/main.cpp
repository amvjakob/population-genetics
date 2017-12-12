#include <gtest/gtest.h>
#include <tclap/CmdLine.h>
#include "../src/Random.hpp"
#include "../src/Data.hpp"
#include "../src/SimulationsExecutor.hpp"

using namespace std;

// For the following tests, we know beforehand the data from the test files given on Moodle
// We test that each data has been well read and collected

TEST(DataReading, PopulationSize) {
	Data data("../data/test_input.txt","../data/test.fa");

	EXPECT_EQ(data.getPopulationSize(), 21);
}

TEST(DataReading, NumberGenerations) {
	Data data("../data/test_input.txt","../data/test.fa");

	EXPECT_EQ(data.getNbGenerations(), 3000);
}

TEST(DataReading, NumberAlleles) {
	Data data("../data/test_input.txt", "../data/test.fa");

	EXPECT_EQ((int) data.getNbAlleles(), 2);
}

TEST(DataReading, MarkerSites) {
	Data data("../data/test_input.txt", "../data/test.fa");

	vector<unsigned int> knownMS = {0, 1, 2, 3};

	EXPECT_EQ(data.getMarkerSites(), knownMS);
}


TEST(DataReading, InitialFrequencies) {
	Data data("../data/test_input.txt","../data/test.fa");

	vector<int> knownCount = {9, 12};

	for (size_t i(0); i < (size_t) data.getAllelesCount().size(); ++i) {
		EXPECT_NEAR(data.getAllelesCount()[i], knownCount[i], 1E-3);
	}
}

TEST(DataReading, NumberReplicates) {
	Data data("../data/test_input.txt","../data/test.fa");

	EXPECT_EQ(data.getNbReplicates(), 500);
}

TEST(DataReading, NucleotidesMutations) {
	Data data("../data/test_input.txt","../data/test.fa");

	vector <double> knownMut = {10E-8, 10E-8, 10E-8, 10E-8};

	for (size_t i(0); i < data.getMutationRates().size() ; ++i) {
		EXPECT_NEAR(data.getMutationRates()[i], knownMut[i], 1E-3);
	}
}

TEST (DataReading, Bottleneck) {
	Data data("../data/test_input.txt","../data/test.fa");
	
	EXPECT_EQ(data.getPopReduction(), 2);
	EXPECT_EQ(data.getBottleneckStart(), 20);
	EXPECT_EQ(data.getBottleneckEnd(), 40);
}

TEST(DataReading, AlleleSelection) {
	Data data("../data/test_input.txt","../data/test.fa");

	vector <double> knownSel = {0.1, -0.8};

	for (size_t i(0); i < data.getSelections().size() ; ++i) {
		EXPECT_NEAR(data.getSelections()[i], knownSel[i], 1E-3);
	}
}

TEST(SelectionTest, AlleleLethality) {
    Data data("../data/test_input.txt","../data/test.fa");

    vector <double> knownProbabilities = {0.5, -1};

    data.collectAll();
    Simulation simul = Simulation({"1", "2"}, {10, 20}, knownProbabilities);

    simul.update(1);

    EXPECT_EQ(simul.getAllelesCount()[1], 0.0);
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

TEST (MigrationTest, FixSubPopulation) {
	
	std::vector<std::string> alleles = { "1", "2", "3" };
	std::vector< std::vector< unsigned int > > subPopulations = { { 10, 0, 0 }, { 0, 20, 0 }, { 0, 0, 30 } };
	std::vector< std::vector< unsigned int > > migrationRates = {
		{ 0, 3, 5 },
		{ 3, 0, 6 },
		{ 5, 6, 0 } 
	};
	
	Simulation simul = Simulation(alleles, subPopulations, migrationRates);
	
	/*
	for (auto& subPop : subPopulations) {
		for (auto& subPopAllele : subPop) {
			std::cout << subPopAllele << '\t';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	
	for (auto& subPop : simul.getSubPopulations()) {
		for (auto& subPopAllele : subPop) {
			std::cout << subPopAllele << '\t';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	
	for (auto& subPop : simul.getSubPopulationSizes()) {
		std::cout << subPop << '\t';
	}
	std::cout << std::endl;
	* 
	* assert(false);
			ASSERT_TRUE(false);
	*/
	

	int t(0);
	while (t < 500) {
		simul.update(t);
		
				
		int totalSize = 0;
		for (size_t i(0); i < simul.getSubPopulations().size(); ++i) {			
			int subPopulationSize = 0;			
			for (auto& subPopAlleleCount : simul.getSubPopulations()[i])
				subPopulationSize += subPopAlleleCount;
			
			// checking if at each step, the incomes are equal to outcomes for each subPop
			EXPECT_EQ(subPopulationSize, simul.getSubPopulationSizes()[i]);
			
			totalSize += subPopulationSize;
		}
		
		EXPECT_EQ(totalSize, simul.getPopulationSize());

		++t;
	}
}

TEST (MigrationTest, MigrationEffect) {
	std::vector<std::string> alleles = { "1", "2", "3" };
	std::vector< std::vector< unsigned int > > subPopulations = { { 10, 0, 0 }, { 0, 20, 0 }, { 0, 0, 30 } };
	std::vector< std::vector< unsigned int > > migrationRates = { 
		{ 0, 3, 5 },
		{ 3, 0, 6 },
		{ 5, 6, 0 } 
	};
	
	Simulation simul = Simulation(alleles, subPopulations, migrationRates);

	int t(0);
	while (t < 500) {
        simul.update(t);
		++t;
	}

	for (size_t i(0); i < simul.getSubPopulations().size(); ++i) {
		for (size_t j(0); j < simul.getSubPopulations()[i].size(); ++j) {
			// making sure that each subpop is different to its initial state --> proof of migrations
            if (i != j) {
                EXPECT_TRUE(subPopulations[i][j] <= simul.getSubPopulations()[i][j]);
            }
		}
	}
}

TEST (MigrationTest, CompleteGraphTest ) {
	
	std::vector<std::string> alleles = { "1", "2", "3" };
	std::vector< std::vector< unsigned int > > subPopulations = { { 10, 0, 0 }, { 0, 20, 0 }, { 0, 0, 30 } };
	std::vector< std::vector< unsigned int > > migrationRates = { 
		{ 0, 3, 5 },
		{ 3, 0, 6 },
		{ 5, 6, 0 } 
	};
	
	Simulation simul = Simulation(alleles, subPopulations, migrationRates);

	int t(0);
	while (t < 500) {
        simul.update(t);
		++t;
	}
	
	/*
    std::cout << "SubPOPTEST" << std::endl;
    for (auto& mig : subPopulations) {
        for (auto& m : mig) {
            std::cout << m << '\t';
        }
        std::cout << std::endl;
    }

    std::cout << "SubPOP" << std::endl;
    for (auto& mig : simul.getSubPopulations()) {
        for (auto& m : mig) {
            std::cout << m << '\t';
        }
        std::cout << std::endl;
    }
    * */
    
    for (size_t i(0); i < simul.getSubPopulations().size(); ++i) {
		for (size_t j(0); j < simul.getSubPopulations()[i].size(); ++j) {
			// making sure that each subpop is different to its initial state --> proof of migrations
            if (i != j) {
                EXPECT_TRUE(subPopulations[i][j] <= simul.getSubPopulations()[i][j]);
            }
		}
	}
}


TEST (MigrationTest, RingTest) {

    std::vector<std::string> alleles = { "1", "2", "3" };
	std::vector< std::vector< unsigned int > > subPopulations = { { 10, 0, 0 }, { 0, 20, 0 }, { 0, 0, 30 } };
	std::vector< std::vector< unsigned int > > migrationRates = { 
		{ 0, 3, 5 },
		{ 3, 0, 6 },
		{ 5, 6, 0 } 
	};
	
	Simulation simul = Simulation(alleles, subPopulations, migrationRates);

	int t(0);
	while (t < 500) {
        simul.update(t);
		++t;
	}
    
    for (size_t i(0); i < simul.getSubPopulations().size(); ++i) {
		for (size_t j(0); j < simul.getSubPopulations()[i].size(); ++j) {
			// making sure that each subpop is different to its initial state --> proof of migrations
            if (i != j) {
                EXPECT_TRUE(subPopulations[i][j] <= simul.getSubPopulations()[i][j]);

            }
		}
	}
}

TEST (MigrationTest, StarTest) {

	std::vector<std::string> alleles = { "1", "2", "3" };
	std::vector< std::vector< unsigned int > > subPopulations = { { 10, 0, 0 }, { 0, 20, 0 }, { 0, 0, 30 } };
	std::vector< std::vector< unsigned int > > migrationRates = { 
		{ 0, 3, 5 },
		{ 3, 0, 6 },
		{ 5, 6, 0 } 
	};
	size_t starCenter = 1;
	
	Simulation simul = Simulation(alleles, subPopulations, migrationRates);

	int t(0);
	while (t < 500) {
        simul.update(t);
		++t;
	}
    
    for (size_t i(0); i < simul.getSubPopulations().size(); ++i) {
		for (size_t j(0); j < simul.getSubPopulations()[i].size(); ++j) {
			// making sure that each subpop is different to its initial state --> proof of migrations
            if (i != j) {
				
				EXPECT_TRUE(subPopulations[i][i] >= simul.getSubPopulations()[i][i]);
				
				if (i == starCenter) {
					// making sure that each subpop exchanges with all the others
					EXPECT_TRUE(subPopulations[i][j] <= simul.getSubPopulations()[i][j]);

				} else {
					EXPECT_TRUE(subPopulations[i][starCenter] <= simul.getSubPopulations()[i][starCenter]);
				}
            }
		}
	}
}

TEST (BottleneckTest, PopulationReduction) {
	int startTime = 20;
	int endTime = 40;
	double reduction = 2.0;
	
	int initialPopulation = 20;
	int reducedPopulation = 10;
	
    Simulation simul = Simulation({ "1", "2" }, { 10, 10 }, startTime, endTime, reduction);
    
    double t(0);
    while(t < 500) {
        simul.update(t);

        if (t < startTime) {
			EXPECT_EQ(simul.getPopulationSize(), initialPopulation);
			
		} else if (t >= startTime and t < endTime) {
			EXPECT_EQ(simul.getPopulationSize(), reducedPopulation);
			
		} else if (t >= endTime) {
			EXPECT_EQ(simul.getPopulationSize(), initialPopulation);
		}
		
		++t;
    }
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
