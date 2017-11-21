#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "SimulationsExecutor.hpp"


SimulationsExecutor::SimulationsExecutor(int n, int populationSize, 
	int simulationSteps, std::vector<double> fqs)
  : isFullMode(false), nSimulations(n), N(populationSize),
	  T(simulationSteps), alleleFqs(fqs)
{
	prepare();
}

SimulationsExecutor::SimulationsExecutor(Data& data)
  : isFullMode(true)
{
	nSimulations = data.getReplicates();
	N = data.getPopSize();
	T = data.getGenerations();
	
	alleleFqs = data.getAlleleFqs();
	
	markerSites = data.getMarkerSites();
	
	sequences = data.getSequences();
	sequences.sort();
	sequences.unique();
	
	mutations = data.getMutations();
	
	// generate allele map
	int i = 0;
	for (auto it = sequences.begin(); it != sequences.end(); ++it) {
		// get initial number of allele in population
		double alleleNb = alleleFqs[i] * N;
		
		// number should be an int, e.g. frequency should be realistic and tied to population size
		assert(alleleNb - ((int) alleleNb) < 1E-4);
		
		std::string idx = *it;
		alleles[idx] = (int) alleleNb;
		
		++i;
	}
	
	prepare();
}


void SimulationsExecutor::prepare() {
	// open result file
	results.open("results.txt");
	
	// init buffer access values
	bufferLowestStep = 0;
	bufferHighestStep = 0;
	
	// add first empty buffer data
	outputBuffer.push_back( std::vector<std::string>(nSimulations) );
}


void SimulationsExecutor::execute() {
	// create correct number of threads
	threads = std::vector<std::thread>(nSimulations);
  
	// init each thread
	for (int i = 0; i < nSimulations; ++i) {
		threads[i] = std::thread(
			[=] { runSimulation(i); }
		);
	}

	// join threads (wait for every thread to end before ending main thread)
	for (auto& th : threads) th.join(); 
}


void SimulationsExecutor::runSimulation(int id) {
	// create new simulation
	Simulation simul = isFullMode ? Simulation(alleles) : Simulation(N, alleleFqs);
	
	std::vector<std::string> states(T + 2);
	
	// write initial allele frequencies
	states[0] = simul.getAlleleFqsForOutput();
	// writeData(simul.getAlleleFqsForOutput(), id, 0);
	
	int t = 0;
	while (t < T) {
		// update simulation
		simul.update();
		
		// increment clock
		++t;
		
		// write allele frequencies
		states[t] = simul.getAlleleFqsForOutput();
		// writeData(simul.getAlleleFqsForOutput(), id, t);
	}
	
	// write final line: allele identifiers
	states[t + 1] = simul.getAlleleStrings();
	// writeData(simul.getAlleleStrings(), id, t + 1);
	
	
	std::size_t lineLength = states.back().size();
	std::size_t precision = simul.getPrecision();
	
	for (auto& state : states) {
		if (state.size() < lineLength) {
			std::stringstream ss;
			ss << state;
			
			while (ss.tellp() < (int) lineLength) {
				ss << '|' << std::setprecision(precision) << std::fixed << 0.0;
			}
		}
	}
	
	for (int i = 0; i < (int) states.size(); ++i) {
		writeData(states[i], id, i);
	}
}

void SimulationsExecutor::writeData(std::string data, int threadId, int step) {
	// check that step is valid
	if (step < bufferLowestStep) { 
		std::string msg = "Error: trying to add data of a step that has already been written to the result file.";
		std::cerr << msg << std::endl;
		throw msg;
	}
	
	// start lock here to restrict the execution of this function to a single thread 
	std::lock_guard<std::mutex> lock(writerMutex);
	
	if (step > bufferHighestStep) {
		// assign new highest step
		bufferHighestStep = step;
		
		// insert new empty vector
		outputBuffer.push_back(std::vector<std::string>(nSimulations));
	}
	
	// add data to buffer
	outputBuffer[step - bufferLowestStep][threadId] = data;
	
	// check for data completeness for the lowest step
	if (std::any_of(
			std::begin(outputBuffer.front()),
			std::end(outputBuffer.front()),
			[](std::string elem) {
				return elem.empty();
			})
		) return;
	
	// write data from buffer to file, since the buffer for the
	// lowest step is full
	writeAlleleFqs(outputBuffer.front());
	outputBuffer.pop_front();
	++bufferLowestStep;
}

void SimulationsExecutor::writeAlleleFqs(const std::vector<std::string>& alleleFqs) {
	for (auto const& data : alleleFqs) {
		results << data << '\t';
	}
	
	results << '\n';
}
