#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
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
		alleles[idx] = (int) alleleNb;
	}
	
	
	int alleleCount = 0;
	for (auto const& allele : alleles)
		alleleCount += allele.second;


	// population size should equal total allele number
	assert(alleleCount == populationSize);
	
	calcOutputConstants();
}

Simulation::Simulation(const std::unordered_map<std::string, int>& als) {
	alleles = als;
	
	int alleleCount = 0;
	for (auto const& allele : alleles)
		alleleCount += allele.second;
		
	populationSize = alleleCount;
	
	calcOutputConstants();
}

void Simulation::calcOutputConstants() {
	std::size_t alleleIdSize = (*(alleles.begin())).first.size();
	
	// the recurring 2 is the size of '0.', the part before the precision
	precision = alleleIdSize <= 2 + _MIN_OUTPUT_PRECISION_ ? _MIN_OUTPUT_PRECISION_ : alleleIdSize - 2;
	additionalSpaces = std::max(precision + 2 - alleleIdSize, (size_t) 0);
}

const std::unordered_map<std::string, int>& Simulation::getAlleles() const {
	return alleles;
}

std::string Simulation::getAlleleFqsForOutput() const {
	std::stringstream ss;
	
	for (auto allele = alleles.begin(); allele != alleles.end(); ++allele) {
		if (allele != alleles.begin()) ss << '|';
		ss << std::setprecision(precision) << std::fixed << (*allele).second * 1.0 / populationSize;
	}
	
	return ss.str();
}

std::string Simulation::getAlleleStrings() const {
	std::stringstream ss;
	
	for (auto allele = alleles.begin(); allele != alleles.end(); ++allele) {
		if (allele != alleles.begin()) ss << '|';
		ss << (*allele).first << std::string(additionalSpaces, ' ');
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
	
	std::cout << "Update" << std::endl;
	
	int alleleCount = 0;
	for (auto const& allele : alleles)
		alleleCount += allele.second;
		
	std::cout << "Pop size: " << alleleCount << std::endl;
	std::cout << "Number of alleles: " << alleles.size() << std::endl;
	
	std::cout << "Alleles" << std::endl;
				for (auto& allele : alleles)
					std::cout << allele.first << ": " << allele.second << std::endl;
	
	for (auto& allele : alleles) {
		// remaining parent population size should be 0 or more	
		assert(nParent >= 0);
		
		std::cout << nParent << std::endl;
		
		double p;
		
		if (nParent == 0) {
			// the only way for the parent population to be 0 is if the 
			// current allele is not present anymore (wiped out)
			assert(allele.second == 0);
			p = 0;
		} else {
			// generate new allele copy number
			p = allele.second * 1.0 / nParent;
		}
		
		// reduce residual "parent population size" / "gene pool"
		nParent -= allele.second;
		
		// generate new number of allele copies in population
        int newAlleleCount = RandomDist::binomial(populationSize - nOffspring, p);
		allele.second = newAlleleCount;
		
		// reduce residual population size
		nOffspring += newAlleleCount;
	}
	
	std::cout << nParent << std::endl;
	
	// after new population generation, residual population size should be 0
	assert(nParent == 0);
	assert(nOffspring == populationSize);
	
	// additional params here!
	switch (executionMode) {
		case _PARAM_MUTATIONS_:
			{
				// todo use specific values for marker sites
				double PROBABILITY = 0.5;
				
					
				std::vector< std::map<Nucleotide, size_t> > markerCount((*(alleles.begin())).first.size());
				std::map<std::string, int> mutatedAlleles;
				
				for (auto& allele : alleles) {
					for (size_t i = 0; i < allele.first.size(); ++i) {
						markerCount[i][Allele::charToNucl.at(allele.first[i])] += allele.second;
					}
				}
				
				for (int markerIdx = 0; markerIdx < (int) markerCount.size(); ++markerIdx) {
					
					for (auto& nuclCount : markerCount[markerIdx]) {
						int mutX = RandomDist::binomial(nuclCount.second, PROBABILITY * 1.0 / nuclCount.second);
						
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
								for (auto& allele : alleles) {
									if (Allele::charToNucl.at(allele.first[markerIdx]) == nuclCount.first) {
										sourceCount += allele.second;
										
										if (source < sourceCount) {
											// we got our source
											
											// create new mutated allele
											std::string newAllele = std::string(allele.first);
											newAllele[markerIdx] = Allele::nuclToChar[(int) target];
											
											// remove original allele
											if (mutatedAlleles.find(allele.first) != mutatedAlleles.end()) {
												mutatedAlleles[allele.first]--;
											} else {
												mutatedAlleles[allele.first] = -1;
											}
											
											// add mutated allele
											if (mutatedAlleles.find(newAllele) != mutatedAlleles.end()) {
												mutatedAlleles[newAllele]++;
											} else {
												mutatedAlleles[newAllele] = 1;
											}
										}
									}								
								}
							}
						}
					}
				}
				
				// merge alleles and mutated alleles
				for (auto mutatedAllele : mutatedAlleles) {
					alleles[mutatedAllele.first];
					alleles[mutatedAllele.first] += mutatedAllele.second;
				}
				
			}
			
			break;
			
		default:
			break;		
	}
}


