#include <gtest/gtest.h>
#include <tclap/CmdLine.h>
#include "../src/Random.hpp"
#include "../src/Data.hpp"
#include "../src/SimulationsExecutor.hpp"

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

TEST(DataReading, AlleleSelection) {
	Data data("../data/test_input.txt","../data/test.fa");

	vector <double> knownSel = {0.1, -0.8};

	data.collectAll();

	for (size_t i(0); i < data.getSelections().size() ; ++i) {
		EXPECT_NEAR(data.getSelections()[i], knownSel[i], 1E-3);
	}
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



TEST (MigrationTest , FixSubPopulation){

	Data data("../data/test_input.txt","../data/test.fa");

	data.collectAll();

	SimulationsExecutor simulationsExecutor(data);

	Simulation simul = simulationsExecutor.createSimulation();

	double stopTime (data.getGenerations());

	double t(0);

	while (t <stopTime ){

		simul.update(t);

		for(size_t i(0); i<simul.getAlleleCount().size();++i  ) {

			//checking if at each step , the incomes are equal to outcomes for each subPop
			EXPECT_EQ(simul.subPopulationSize(simul.getSubPop()[i]),simul.getAlleleCount()[i]);
		}

		++t;
	}
}

TEST (MigrationTest ,MigrationEffect) {

	Data data("../data/test_input.txt","../data/test.fa");

	data.collectAll();

	SimulationsExecutor simulationsExecutor(data);

	Simulation simul = simulationsExecutor.createSimulation();

	double stopTime (data.getGenerations());

	double t(0);

	//to initialize the comparison vector
	std::vector< std::vector<unsigned int> > subPopTest = simul.getSubPop();


	while (t <stopTime ){

        simul.update((int)t);
		++t;
	}

	for(size_t i(0); i<simul.getSubPop().size();++i  ) {


		for(size_t j(0); j<simul.getSubPop()[i].size();++j) {
			//making sure that each subpop is different to it initial state --> proof of migrations
            if (i!=j) {
                EXPECT_TRUE(subPopTest[i][j] <= simul.getSubPop()[i][j]);
            }
		}

	}

}

TEST (MigrationTest , CompleteGrapshTest ) {


    Data data("../data/test_input.txt", "../data/test.fa");

    data.collectAll();

    data.setDataMigTest(_PARAM_MIGRATION_,_COMPLETE_GRAPH_,_RANDOM_);

    SimulationsExecutor simulationsExecutor(data);

    Simulation simul = simulationsExecutor.createSimulation();

    double stopTime(data.getGenerations());

    double t(0);

    //to initialize the comparison subPopulation  vector

	std::vector<std::vector<unsigned int> > subPopTest = simul.getSubPop();

    while (t < stopTime) {

        simul.update(t);

        ++t;

    }

    std::cout << "SubPOPTEST" << std::endl;
    for (auto& mig : subPopTest) {
        for (auto& m : mig) {
            std::cout << m << '\t';
        }
        std::cout << std::endl;
    }

    std::cout << "SubPOP" << std::endl;
    for (auto& mig : simul.getSubPop()) {
        for (auto& m : mig) {
            std::cout << m << '\t';
        }
        std::cout << std::endl;
    }

    for (size_t i(0); i < simul.getSubPop().size(); ++i) {


        for (size_t j(0); j < simul.getSubPop()[i].size(); ++j) {


            if (i!=j) {

                //making sure that each subpop exchanges with all the others

                EXPECT_TRUE(subPopTest[i][i] >= simul.getSubPop()[i][i]);

                EXPECT_TRUE(subPopTest[i][j] <= simul.getSubPop()[i][j]);

            }
        }
    }
}


TEST (MigrationTest , RingTest ) {


    Data data("../data/test_input.txt", "../data/test.fa");

    data.collectAll();

    data.setDataMigTest(_PARAM_MIGRATION_,_RING_,_RANDOM_);

    SimulationsExecutor simulationsExecutor(data);

    Simulation simul = simulationsExecutor.createSimulation();


    double stopTime(data.getGenerations());

    double t(0);

    //to initialize the comparison subPopulation  vector
    std::vector<std::vector<unsigned int> > subPopTest = simul.getSubPop();

    while (t < stopTime) {
        simul.update(t);
        ++t;
    }

    std::cout << "SubPOPTEST" << std::endl;
    for (auto& mig : subPopTest) {
        for (auto& m : mig) {
            std::cout << m << '\t';
        }
        std::cout << std::endl;
    }

    std::cout << "SubPOP" << std::endl;
    for (auto& mig : simul.getSubPop()) {
        for (auto& m : mig) {
            std::cout << m << '\t';
        }
        std::cout << std::endl;
    }


    for (size_t i(0); i < simul.getSubPop().size(); ++i) {


        for (size_t j(0); j < simul.getSubPop()[i].size(); ++j) {

            if (j == i + 1) {

                //making sure that each subpop exchanges with all the others

                EXPECT_TRUE(subPopTest[i][i] >= simul.getSubPop()[i][i]);

                EXPECT_TRUE(subPopTest[i][j] <= simul.getSubPop()[i][j]);

            }

        }
    }
}

TEST (MigrationTest , StarTest ) {


    Data data("../data/test_input.txt", "../data/test.fa");

    data.collectAll();

    data.setDataMigTest(_PARAM_MIGRATION_,_STAR_,_RANDOM_);


    SimulationsExecutor simulationsExecutor(data);

    Simulation simul = simulationsExecutor.createSimulation();

    double stopTime(data.getGenerations());

    double t(0);

    //to initialize the comparison subPopulation  vector
    std::vector<std::vector<unsigned int> > subPopTest = simul.getSubPop();

    while (t < stopTime) {

        simul.update(t);

        ++t;

    }

    std::cout << "SubPOPTEST" << std::endl;
    for (auto& mig : subPopTest) {
        for (auto& m : mig) {
            std::cout << m << '\t';
        }
        std::cout << std::endl;
    }

    std::cout << "SubPOP" << std::endl;
    for (auto& mig : simul.getSubPop()) {
        for (auto& m : mig) {
            std::cout << m << '\t';
        }
        std::cout << std::endl;
    }

    for (size_t i(0); i < simul.getSubPop().size(); ++i) {


        for (size_t j(0); j < simul.getSubPop()[i].size(); ++j) {

            if (i != j) {

						if (i == simulationsExecutor.getStarCenter()) {

							//making sure that each subpop exchanges with all the others

							EXPECT_TRUE(subPopTest[i][i] >= simul.getSubPop()[i][i]);

							EXPECT_TRUE(subPopTest[i][j] <= simul.getSubPop()[i][j]);

						} else {

							EXPECT_TRUE(subPopTest[i][i] >= simul.getSubPop()[i][i]);

							EXPECT_TRUE(subPopTest[i][simulationsExecutor.getStarCenter()] <=
										simul.getSubPop()[i][simulationsExecutor.getStarCenter()]);
						}
            }

        }

    }


}


TEST (BottleneckTest, PopulationReduction) {
	
	Data data("../data/test_input.txt", "../data/test.fa");
	
    data.collectAll();
	
	data.setExecutionMode(_PARAM_BOTTLENECK_);
	
    SimulationsExecutor simulationsExecutor(data);
	
    Simulation simul = simulationsExecutor.createSimulation();
    
    double stopTime(data.getGenerations());
	
	double reduction = data.getPopReduction();
	
	int population = data.getPopSize();
	
	int popReducted = population / reduction;
	
    double t(1);

	
    while (t < stopTime) {
		
        simul.update(t);

        ++t;

        if (t <= data.getBottleneckStart()) {
			EXPECT_EQ(simul.getPopSize(), population);	
		} else if (t>=data.getBottleneckStart() and t <=data.getBottleneckEnd()) {
			EXPECT_EQ(simul.getPopSize(),popReducted);
		} else if (t > data.getBottleneckEnd()) {
			if (population%2==0) {
			EXPECT_EQ(simul.getPopSize(), population);
			} else {
			EXPECT_EQ(simul.getPopSize(), population-1);
			}
		}
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
