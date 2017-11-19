#ifndef SIMULATIONS_EXECUTOR_H
#define SIMULATIONS_EXECUTOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <mutex>
#include <thread>
#include "Simulation.hpp"

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
	 * \param alleleFqs				A vector of initial allele frequenciess
	 * */
	SimulationsExecutor(int n, int populationSize, int simulationSteps, std::vector<double> alleleFqs);
	
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
	 * \param alleleFqs		a vector of strings (each being a formatted list of allele frequencies)
	 * */
	void writeAlleleFqs(const std::vector<std::string>& alleleFqs);
	

private:
	
	//!< Number of simulations to be executed simultaneously
	int nSimulations;
	
	//!< Size of the population in a Simulation
	int N;
	
	//!< Length of a Simulation in steps
	int T;
	
	//!< Initial allele frequencies for a Simulation
	std::vector<double> alleleFqs;
	
	
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
};

#endif
