#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip> 
#include "Simulation.hpp"

Simulation::Simulation(int N, int T, std::vector<double> alleleFq)
	: populationSize(N), simulationSteps(T)
{
	// make sure sensible parameters were used
	assert(N > 0);
	assert(T > 0);
	
	for (int i = 0; i < (int) alleleFq.size(); ++i) {
		// get initial number of allele in population
		double alleleNb = alleleFq[i] * populationSize;
		
		// number should be an int, e.g. frequency should be realistic and tied to population size
		assert(alleleNb - ((int) alleleNb) < 1E-4);
		
		alleles[i] = (int) alleleNb;
	}
	
	
	int alleleCount = 0;
	for (auto const& allele : alleles)
		alleleCount += allele.second;


	// population size should equal total allele number
	assert(alleleCount == populationSize);
}

const std::map<int, int>& Simulation::getAlleles() const {
	return alleles;
}

std::string Simulation::getAlleleFqsForOutput() const {
	std::stringstream ss;
	
	for (auto allele = alleles.begin(); allele != alleles.end(); ++allele) {
		ss << std::setprecision(2) << std::fixed << (*allele).second * 1.0 / populationSize;
		if (allele != alleles.end()) ss << '|';
	}
	
	return ss.str();
}

std::string Simulation::getAlleleStrings() const {
	std::stringstream ss;
	
	for (auto const& allele : alleles) {
		ss << allele.first << "   |";
	}
	
	return ss.str();
}

void Simulation::update(RandomDist& randomDist) {
	int nParent = populationSize;
	int nOffspring = 0;
	
	for (auto& allele : alleles) {
		// remaining parent population size should be 0 or more	
		assert(nParent >= 0);
		
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
        int newAlleleCount = randomDist.binomial(populationSize - nOffspring, p);
		allele.second = newAlleleCount;
		
		// reduce residual population size
		nOffspring += newAlleleCount;
	}
	
	// after new population generation, residual population size should be 0
	assert(nParent == 0);
	assert(nOffspring == populationSize);
}


