#ifndef SIMULATIONS_EXECUTOR_H
#define SIMULATIONS_EXECUTOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <map>
#include <unordered_map>
#include <mutex>
#include <thread>
#include "Simulation.hpp"
#include "Data.hpp"
#include "Allele.hpp"

/** \brief Class representing a SimulationsExecutor
 * 
 * This class is a wrapper for the execution of multiple single
 * Simulations with the same initial parameters (e.g. to generate statistics).
 * The Simulations are run on different threads and their outputs stored
 * in a common file. * 
 * */
class SimulationsExecutor {
	
public:

	/** \brief SimulationsExecutor constructor
	 * 
	 * Initialises a new series of Simulations
	 * 
	 * \param n						n, the number of "identical" simulations to be run
	 * \param populationSize		N, the number of individuals in the population
	 * \param simulationSteps		T, the number of generations in the simulation
	 * \param allelesCount			A vector of initial allele counts
	 * */
	SimulationsExecutor(int n, int populationSize, int simulationSteps, std::vector<unsigned int> allelesCount);
	
	/** \brief SimulationsExecutor constructor
	 * 
	 * Initialises a new series of Simulations with a given Data object
	 * 
	 * \param data					Data file
	 * */
	SimulationsExecutor(const Data& data);
	
	
	//!< Because of the threads, we do not allow any copies
	SimulationsExecutor(const SimulationsExecutor& other) = delete;
	
	//!< Because of the threads, we do not allow any copies
	SimulationsExecutor& operator=(const SimulationsExecutor& other) = delete;
	
	/** \brief Start the execution of the simulations
	 * 
	 * Invoking this method will create the threads and execute them.
	 * */
	void execute();
	
	
protected: 

	/** \brief Run a simulation
	 * 
	 * This method is executed by a thread.
	 * 
	 * \param id		the id of the simulation
	 * */
	void runSimulation(int id);
	
	/** \brief Write data to a buffer to be eventually put in a result file
	 * 
	 * Writes data to a temporary buffer, which performs the actual output
	 * once it has gathered the data of every thread for a given step.
	 * This function is lock guarded.
	 * 
	 * \param data			the data to be written (usually allele frequencies)
	 * \param threadId		the id of the executing thread
	 * \param step			the current step of the simulation
	 * */
	void writeData(std::string data, int threadId, int step);
	
	/** \brief Write one step of all simulations to the result file
	 * 
	 * Wrties the data passed as argument in the result file.
	 * 
	 * \param step			the step number of the simulation to be written
	 * \param alleleFqs		a vector of strings (each being a formatted list of allele frequencies)
	 * */
	void writeAlleleFqs(int step, const std::vector<std::string>& alleleFqs);
    
    
    void writeDataMigs(std::string data, int threadId, int step);

	
	/** \brief Generate table for nucleotide mutation rates based on user input data
	 * 
	 * Automatically selects the correct model based on user data
	 * */
	 void generateMutationRates(const Data& data);
	 
	 /** \brief Generate the subpopulations for a Simulation with migration
	  * */
	 void generateSubPopulations(const Data& data);
    
    
    void writeAlleleMigFqs(const std::vector<std::string>& alleleFqs);

	

private:

	//!< Prepare SimulationExecutor
	void prepare();
	
	
	//!< Execution mode of simulation
	const int executionMode;
	
	//!< Number of simulations to be executed simultaneously
	int nSimulations;
	
	//!< Size of the population in a Simulation
	int N;
	
	//!< Length of a Simulation in steps
	int T;
	
	//!< List of alleles int the Simulations
	std::vector<std::string> alleles;
	
	//!< Initial allele counts for a Simulation
	std::vector<unsigned int> allelesCount;

	//!< Vector of double containing the user marker sites
	std::vector<double> markerSites;

	//!< Vector of double containing the mutations probabilities of the marker sites
	std::vector<double> mutations;
	
	//!< Table of mutation probabilities
	std::array< std::array<double, Nucleotide::N >, Nucleotide::N > nuclMutationProbs;
	
	
	//!< List of selections frequencies of each alleles
	std::vector<double> selectionFqs;
	
	
	//!< Table containing the sub-populations
    std::vector< std::vector<unsigned int> > subPopulations;
    
    //!< Table containing migration rates for each sub group
    std::vector< std::vector<unsigned int> > migrationRates;
	
	
	
	//!< Result file
	std::ofstream results;
   
    
	//!< Mutex for lock guarding
	std::mutex writerMutex;
	
	//!< List of all threads to be executed
	std::vector<std::thread> threads;
	
	//!< Output buffer for result data
	std::deque< std::vector<std::string> > outputBuffer;

	//!< Lowest step that is still in the buffer
	int bufferLowestStep;
	//!< Highest step that is already in the buffer
	int bufferHighestStep;
	
	
    std::mutex writerMutex2;
	
    std::ofstream migrationResults;
    
    std::deque< std::vector<std::string> > outputBufferMigs;
    //!< Lowest step that is still in the buffer
    int buffer2LowestStep;
    //!< Highest step that is already in the buffer
    int buffer2HighestStep;
};

#endif
