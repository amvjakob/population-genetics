#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "SimulationsExecutor.hpp"


SimulationsExecutor::SimulationsExecutor(int n, int populationSize, 
	int simulationSteps, std::vector<double> fqs)
  : isFullMode(false), nSimulations(n), N(populationSize),
	  T(simulationSteps), alleleFqs(fqs), executionMode(_PARAM_NONE_)
{
	prepare();
}

SimulationsExecutor::SimulationsExecutor(const Data& data)
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
	
	executionMode = data.getExecutionMode();
	
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
	
	// generate nucleotide mutation probabilities according to model
	if (executionMode == _PARAM_MUTATIONS_) {
		std::cout << "Running with param mutations" << std::endl;
		generateMutationRates(data);
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
	Simulation simul = isFullMode ? Simulation(alleles, executionMode, mutations, nuclMutationProbs) : Simulation(N, alleleFqs);
	
	std::vector<std::string> states(T + 2);
	
	// write initial allele frequencies
	states[0] = simul.getAlleleFqsForOutput();
	
	int t = 0;
	while (t < T) {
		// update simulation
		simul.update();
		
		//bottleneck effect
		simul.bottleneck(t);
		
		// increment clock
		++t;
		
		// write allele frequencies
		states[t] = simul.getAlleleFqsForOutput();
	}
	
	// write final line: allele identifiers
	states[t + 1] = simul.getAlleleStrings();
	
	
	std::size_t lineLength = states.back().size();
	std::size_t precision = simul.getPrecision();
	
	for (auto& state : states) {
		if (state.size() < lineLength) {
			std::stringstream ss;
			ss << state;
			
			while (ss.tellp() < (int) lineLength) {
				ss << '|' << std::setprecision(precision) << std::fixed << 0.0;
			}
			
			state = ss.str();
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
	writeAlleleFqs(bufferLowestStep, outputBuffer.front());
	outputBuffer.pop_front();
	++bufferLowestStep;
}

void SimulationsExecutor::writeAlleleFqs(int step, const std::vector<std::string>& alleleFqs) {
	results << step << '\t';
	
	for (auto const& data : alleleFqs) {
		results << data << '\t';
	}
	
	results << '\n';
}

void SimulationsExecutor::generateMutationRates(const Data& data) {
	switch (data.getMutationModel()) {
		case _MUTATION_MODEL_CANTOR_:
			{
				std::cout << "Cantor model" << std::endl;
				
				double p = 1.0 / 3.0;
			
				nuclMutationProbs = { {
					{ { 0.0, p, p, p } },
					{ { p, 0.0, p, p } },
					{ { p, p, 0.0, p } },
					{ { p, p, p, 0.0 } }
				} };
			}
			break;
			
		case _MUTATION_MODEL_KIMURA_:
			{
				std::cout << "Kimura model" << std::endl;
			
				double transition = data.getKimuraDelta();
				double transversion = (1.0 - transition) / 2.0;
				
				nuclMutationProbs = { {
					//  A    C			   G		   T
					{ { 0.0, transversion, transition, transversion } },
					{ { transversion, 0.0, transversion, transition } },
					{ { transition, transversion, 0.0, transversion } },
					{ { transversion, transition, transversion, 0.0 } }
				} };
			}
			break;
			
		case _MUTATION_MODEL_FELSENSTEIN_:
			{
				std::cout << "Felsenstein model" << std::endl;
			
				std::vector<double> consts = data.getFelsensteinConstants();
				
				for (auto& c : consts) {
					assert(c != 1.0);
				}
				
				double pA = consts[Nucleotide::A] / (1.0 - consts[Nucleotide::A]);
				double pC = consts[Nucleotide::C] / (1.0 - consts[Nucleotide::C]);
				double pG = consts[Nucleotide::G] / (1.0 - consts[Nucleotide::G]);
				double pT = consts[Nucleotide::T] / (1.0 - consts[Nucleotide::T]);
				
				nuclMutationProbs = { {
					{ { 0.0, pC, pG, pT } },
					{ { pA, 0.0, pG, pT } },
					{ { pA, pC, 0.0, pT } },
					{ { pA, pC, pG, 0.0 } }
				} };
			}
			break;
			
		default:
			std::cout << "No model" << std::endl;
			break;
	}
	
	// cout the prob matrix
	for (auto& pX : nuclMutationProbs) {
		for (auto& p : pX) {
			std::cout << p << '\t';
		}
		std::cout << std::endl;
	}
}
