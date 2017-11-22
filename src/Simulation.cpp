#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "Simulation.hpp"
#include "Random.hpp"


std::array< std::array<double, N >, N > Simulation::mutationTable = { {
	{ {0.0, 1.0/3.0, 1.0/3.0, 1.0/3.0} },
	{ {1.0/3.0, 0.0, 1.0/3.0, 1.0/3.0} },
	{ {1.0/3.0, 1.0/3.0, 0.0, 1.0/3.0} },
	{ {1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0} }
} };

Simulation::Simulation(int N, std::vector<double> alleleFq)
	: populationSize(N)
{
	// make sure sensible parameters were used
	assert(N > 0);
	
	for (int i = 0; i < (int) alleleFq.size(); ++i) {
		// get initial number of allele in population
		double alleleNb = alleleFq[i] * populationSize;
		
		// number should be an int, e.g. frequency should be realistic and tied to population size
		assert(alleleNb - ((int) alleleNb) < 1E-4);
		
		std::string idx = std::to_string(i);
		
		alleles.push_back(idx);
		allelesCount.push_back((int) alleleNb);
	}
	
	
	int alleleCount = 0;
	for (auto const& allele : allelesCount)
		alleleCount += allele;


	// population size should equal total allele number
	assert(alleleCount == populationSize);
	
	calcOutputConstants();
}

Simulation::Simulation(const std::unordered_map<std::string, int>& als, const std::vector<double>& mutationRates)
  : mutationFqs(mutationRates)
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
		if (allele != allelesCount.begin()) ss << '|';
		ss << std::setprecision(precision) << std::fixed << (*allele) * 1.0 / populationSize;
	}
	
	return ss.str();
}

std::string Simulation::getAlleleStrings() const {
	std::stringstream ss;
	
	for (auto allele = alleles.begin(); allele != alleles.end(); ++allele) {
		if (allele != alleles.begin()) ss << '|';
		ss << (*allele) << std::string(additionalSpaces, ' ');
	}
	
	return ss.str();
}

std::size_t Simulation::getPrecision() const {
	return precision;
}

void Simulation::update() {
	// genetic drift
	int nParent = populationSize;
	int nOffspring = 0;
	
	assert(alleles.size() == allelesCount.size());
	
	for (auto& allele : allelesCount) {
		// remaining parent population size should be 0 or more	
		assert(nParent >= 0);
		
		double p;
		
		if (nParent == 0) {
			// the only way for the parent population to be 0 is if the 
			// current allele is not present anymore (wiped out)
			assert(allele == 0);
			p = 0;
		} else {
			// generate new allele copy number
			p = allele * 1.0 / nParent;
		}
		
		// reduce residual "gene pool"
		nParent -= allele;
		
		// generate new number of allele copies in population
        int newAlleleCount = RandomDist::binomial(populationSize - nOffspring, p);
		allele = newAlleleCount;
		
		// reduce residual population size
		nOffspring += newAlleleCount;
	}
	
	// after new population generation, residual population size should be 0
	assert(nParent == 0);
	assert(nOffspring == populationSize);
	
	// additional params here!
	switch (executionMode) {
		case _PARAM_MUTATIONS_:
			{							
				std::size_t nbMarkers = alleles.front().size();
				
				for (std::size_t markerIdx = 0; markerIdx < nbMarkers; ++markerIdx) {
					
					std::map<Nucleotide, std::size_t> markerCount;
					
					for (std::size_t i= 0; i < alleles.size(); ++i) {
						markerCount[Allele::charToNucl.at(alleles[i][markerIdx])] += allelesCount[i];
					}
					
					for (auto& nuclCount : markerCount) {						
						int mutX = nuclCount.second > 0 ? RandomDist::binomial(nuclCount.second, mutationFqs[markerIdx] * 1.0 / nuclCount.second) : 0;
						
						if (mutX > 0) {
							
							for (int mutation = 0; mutation < mutX; ++mutation) {
								// right now, simplest model with p = 1/3 for each transition
								double mut = RandomDist::uniformDoubleSingle(0, 1);
								Nucleotide target = N;
								
								double pCount = 0.0;
								for (int i = 0; i < (int) N; ++i) {
									pCount += Simulation::mutationTable[(int) nuclCount.first][i];
									if (mut <= pCount) {
										target = (Nucleotide) i;
										break;
									}							
								}
								
								// find source
								int source = RandomDist::uniformIntSingle(0, nuclCount.second);
								
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
													newAllele)
												);
											 
											if (newAlleleIdx < alleles.size()) {
												allelesCount[newAlleleIdx]++;
											} else {
												alleles.push_back(newAllele);
												allelesCount.push_back(1);
											}
											
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
			
			break;
			
		default:
			break;		
	}
}


