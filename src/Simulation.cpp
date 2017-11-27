#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "Simulation.hpp"
#include "Random.hpp"

Simulation::Simulation(int N, std::vector<int> counts)
	: populationSize(N), executionMode(_PARAM_NONE_)
{
	// make sure sensible parameters were used
	assert(populationSize > 0);
	
	for (std::size_t i = 0; i < counts.size(); ++i) {
		std::string idx = std::to_string(i);
		
		alleles.push_back(idx);
		allelesCount.push_back(counts[i]);
	}
	
	int alleleCount = 0;
	for (auto const& allele : allelesCount)
		alleleCount += allele;

	// population size should equal total allele number
	assert(alleleCount == populationSize);
	
	calcOutputConstants();
}

Simulation::Simulation(const std::unordered_map<std::string, int>& als,
			const int execMode,
			const std::vector<double>& mutationRates, const std::array< std::array<double, N>, N >& nuclMutationProbs)
  : executionMode(execMode), mutationFqs(mutationRates), mutationTable(nuclMutationProbs)
{
	for (auto const& allele : als) {
		alleles.push_back(allele.first);
		allelesCount.push_back(allele.second);
	}
	
	int alleleCount = 0;
	for (auto const& allele : allelesCount)
		alleleCount += allele;
		
	populationSize = alleleCount;
	
	// mutation rates - sanitize input
	while (mutationFqs.size() < alleles.front().size()) {
		mutationFqs.push_back(_DEFAULT_MUTATION_RATE_);
	}
	
	calcOutputConstants();
}

void Simulation::calcOutputConstants() {
	std::size_t alleleIdSize = alleles.front().size();
	
	// the recurring 2 is the size of '0.', the part before the precision
	precision = alleleIdSize <= 2 + _MIN_OUTPUT_PRECISION_ ? _MIN_OUTPUT_PRECISION_ : alleleIdSize - 2;
	additionalSpaces = std::max(precision + 2 - alleleIdSize, (size_t) 0);
}

const std::vector<std::string>& Simulation::getAlleles() const {
	return alleles;
}

const std::vector<unsigned int>& Simulation::getAllelesCount() const {
	return allelesCount;
}

std::string Simulation::getAlleleFqsForOutput() const {
	std::stringstream ss;
	
	for (auto allele = allelesCount.begin(); allele != allelesCount.end(); ++allele) {
		if (allele != allelesCount.begin()) ss << _OUTPUT_SEPARATOR_;
		ss << std::setprecision(precision) << std::fixed << (*allele) * 1.0 / populationSize;
	}
	
	return ss.str();
}

std::string Simulation::getAlleleStrings() const {
	std::stringstream ss;
	
	for (auto allele = alleles.begin(); allele != alleles.end(); ++allele) {
		if (allele != alleles.begin()) ss << _OUTPUT_SEPARATOR_;
		ss << (*allele) << std::string(additionalSpaces, ' ');
	}
	
	return ss.str();
}

std::size_t Simulation::getPrecision() const {
	return precision;
}

void Simulation::update() {
	// genetic drift
	assert(alleles.size() == allelesCount.size());
	
	switch (executionMode) {
		case _PARAM_MUTATIONS_:
			RandomDist::multinomial(allelesCount, populationSize);
			
			mutatePopulation();
			
			break;
			
		case _PARAM_NONE_:
		default:
			RandomDist::multinomial(allelesCount, populationSize);
			
	}
}

void Simulation::mutatePopulation() {
	assert(!mutationFqs.empty());
		
	// get number of marker sites					
	std::size_t nbMarkers = alleles.front().size();
	
	// iterate over all marker sites
	for (std::size_t markerIdx = 0; markerIdx < nbMarkers; ++markerIdx) {
		
		// count the number of each nucleotide for the current marker site
		std::map<Nucleotide, std::size_t> markerCount;
		for (std::size_t i= 0; i < alleles.size(); ++i) {
			markerCount[Allele::charToNucl.at(alleles[i][markerIdx])] += allelesCount[i];
		}
		
		// iterate over all nucleotide counts of the current marker site
		for (auto& nuclCount : markerCount) {	
			
			// generate number of mutations for nucleotide X				
			int mutX = nuclCount.second > 0 ? RandomDist::binomial(nuclCount.second, mutationFqs[markerIdx] * 1.0 / nuclCount.second) : 0;
			
			if (mutX > 0) {				
				// iterate over all mutations
				for (int mutation = 0; mutation < mutX; ++mutation) {
											
					// generate the target mutation
					double mut = RandomDist::uniformDoubleSingle(0.0, 1.0);
					Nucleotide target = N;
					
					double pCount = 0.0;
					for (int i = 0; i < (int) Nucleotide::N; ++i) {
						pCount += mutationTable[(int) nuclCount.first][i];
						if (mut <= pCount) {
							target = (Nucleotide) i;
							break;
						}							
					}
					
					if (target == Nucleotide::N) {
						std::cerr << _ERROR_MUTATION_TARGET_UNFINDABLE_ << std::endl;
						throw _ERROR_MUTATION_TARGET_UNFINDABLE_;
					}
				
					// find source
					int source = RandomDist::uniformIntSingle(0, nuclCount.second - 1);
					
					int sourceCount = 0.0;
					for (std::size_t i = 0; i < alleles.size(); ++i) {
						if (Allele::charToNucl.at(alleles[i][markerIdx]) == nuclCount.first) {
							sourceCount += allelesCount[i];
							
							if (source < sourceCount) {
								// we got our source
								
								// create new mutated allele
								std::string newAllele = std::string(alleles[i]);
								newAllele[markerIdx] = Allele::nuclToChar[(int) target];
								
								// remove original allele
								allelesCount[i]--;
								
								// add mutated allele
								std::size_t newAlleleIdx = std::distance(
									alleles.begin(), 
									std::find(
										alleles.begin(),
										alleles.end(),
										newAllele
									));
								 
								if (newAlleleIdx < alleles.size()) {
									allelesCount[newAlleleIdx]++;
								} else {
									alleles.push_back(newAllele);
									allelesCount.push_back(1);
								}
								
								// some user info
								std::cout << "Mutated a " << alleles[i] << " to a " << newAllele << std::endl;
								
								break;
							}
						}								
					}
				}
			}
		}
	}
}


