#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "Simulation.hpp"
#include "Random.hpp"

Simulation::Simulation(const Simulation& other)
  : executionMode(other.executionMode),
	populationSize(other.populationSize),
	alleles(other.alleles),
	allelesCount(other.allelesCount),
	mutationFqs(other.mutationFqs),
	mutationTable(other.mutationTable),
	selectionFqs(other.selectionFqs),
	subPopulations(other.subPopulations),
	subPopulationSizes(other.subPopulationSizes),
	migrationRates(other.migrationRates),
	precision(other.precision),
	additionalSpaces(other.additionalSpaces)
{	}
	
Simulation& Simulation::operator=(const Simulation& other) {
	executionMode = other.executionMode;
	populationSize = other.populationSize;
	alleles = other.alleles;
	allelesCount = other.allelesCount;
	mutationFqs = other.mutationFqs;
	mutationTable = other.mutationTable;
	selectionFqs = other.selectionFqs;
	subPopulations = other.subPopulations;
	subPopulationSizes = other.subPopulationSizes;
	migrationRates = other.migrationRates;
	precision = other.precision;
	additionalSpaces = other.additionalSpaces;
	
	return *this;
}


Simulation::Simulation(std::vector<unsigned int> counts)
  : executionMode(_PARAM_NONE_), populationSize(0)
{	
	for (std::size_t i = 0; i < counts.size(); ++i) {
		std::string idx = std::to_string(i);
		
		alleles.push_back(idx);
		allelesCount.push_back(counts[i]);
		
		populationSize += counts[i];
	}

	// make sure sensible parameters were used
	assert(populationSize > 0);
	
	calcOutputConstants();
}

Simulation::Simulation(const std::vector<std::string>& als,
						const std::vector<unsigned int>& alsCount,
						const std::vector<double>& mutationRates, 
						const std::array< std::array<double, Nucl::Nucleotide::N>, Nucl::Nucleotide::N >& nuclMutationProbs)
  : executionMode(_PARAM_MUTATIONS_),
	populationSize(0), 
	alleles(als), allelesCount(alsCount), 
	mutationFqs(mutationRates), mutationTable(nuclMutationProbs)
{	
	assert(alleles.size() == allelesCount.size());
	
	for (auto& count : allelesCount)
		populationSize += count;
		
	// make sure sensible parameters were used
	assert(populationSize > 0);
	
	// mutation rates - sanitize input
	while (mutationFqs.size() < alleles.front().size()) {
		mutationFqs.push_back(_DEFAULT_MUTATION_RATE_);
	}
	
	calcOutputConstants();
}

Simulation::Simulation(const std::vector<std::string>& als,
						const std::vector<unsigned int>& alsCount,
						const std::vector<double>& selectionRates)
  : executionMode(_PARAM_SELECTION_),
	populationSize(0),
	alleles(als), allelesCount(alsCount), 
	selectionFqs(selectionRates)
{
	assert(alleles.size() == allelesCount.size());
	
	for (auto& count : allelesCount)
		populationSize += count;
		
	// make sure sensible parameters were used
	assert(populationSize > 0);
	
	// selection rates - sanitize input
	for (auto& sfq : selectionFqs) {
		assert(sfq >= -1.0);
	}
	while (selectionFqs.size() < alleles.size()) {
		selectionFqs.push_back(0.0);
	}
	
	calcOutputConstants();
}
    
Simulation::Simulation(const std::vector<std::string>& als,
						const std::vector< std::vector<unsigned int> >& subPopsCount,
						const std::vector< std::vector<unsigned int> >& migrationFqs)
  : executionMode(_PARAM_MIGRATION_), 
	populationSize(0), 
	alleles(als), subPopulations(subPopsCount),
	migrationRates(migrationFqs)
{    
    for (auto& population : subPopsCount) {
		assert(population.size() == alleles.size());
		
		unsigned int subPopSize = 0;
		for (auto& count : population)
			subPopSize += count;
			
		subPopulationSizes.push_back(subPopSize);
		populationSize += subPopSize;
	}
		
	// make sure sensible parameters were used
	assert(populationSize > 0);
	assert(migrationRates.size() == subPopulations.size());
	assert(migrationRates.front().size() == subPopulations.size());
	
	/*for (size_t i = 0; i < subPopulations.size(); ++i) {
		int outgoing = 0;
		for (auto& out : migrationRates[i])
			outgoing += out;
			
		assert(outgoing <= (int) subPopulationSizes[i]);
	}*/
	
	calcOutputConstants();
}

Simulation::Simulation(const std::vector<std::string>& als,
						const std::vector<unsigned int>& alsCount,
						const int start,
						const int stop,
						const double reduction)
  : executionMode(_PARAM_BOTTLENECK_),
	populationSize(0), 
	alleles(als), allelesCount(alsCount),
	bottleneckStart(start),
	bottleneckEnd(stop),
	popReduction(reduction)
{
	assert(alleles.size() == allelesCount.size());
	
	for (auto& count : allelesCount)
		populationSize += count;
		
	// make sure sensible parameters were used
	assert(populationSize > 0);
	
	calcOutputConstants();
	
	assert(popReduction!=0);
	assert(bottleneckStart<=bottleneckEnd);
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
	
	if (executionMode != _PARAM_MIGRATION_) {
		
		for (auto allele = allelesCount.begin(); allele != allelesCount.end(); ++allele) {
			if (allele != allelesCount.begin()) ss << _OUTPUT_SEPARATOR_;
			ss << std::setprecision(precision) << std::fixed << (*allele) * 1.0 / populationSize;
		}
		
	} else {
		
		if (_MIGRATION_DETAILED_OUTPUT_) {

			for (auto subPop = subPopulations.begin(); subPop != subPopulations.end(); ++subPop) {
				for (auto allele = subPop->begin(); allele != subPop->end(); ++allele) {
					if (allele != subPop->begin()) ss << _OUTPUT_SEPARATOR_;
					ss << std::setprecision(precision) << std::fixed << (*allele) * 1.0 / populationSize;
				}
				
				ss << _MIGRATION_OUTPUT_SEPARATOR_;
			}
		} else {


			int nAlleles = subPopulations.front().size();
			for (int i = 0; i < nAlleles; ++i) {
				int sum = 0;
				for (int j = 0; j < (int) subPopulations.size(); ++j) {
					sum += subPopulations[j][i];
				}
				
				if (i != 0) ss << _OUTPUT_SEPARATOR_;
				ss << std::setprecision(precision) << std::fixed << sum * 1.0 / populationSize;
			}
		}
	}
	
	return ss.str();
}

std::string Simulation::getAlleleStrings() const {
	std::stringstream ss;

	for (auto allele = alleles.begin(); allele != alleles.end(); ++allele) {
		if (allele != alleles.begin()) ss << _OUTPUT_SEPARATOR_;
		ss << (*allele) << std::string(additionalSpaces, ' ');
	}
	
	// add string identifiers for each subpopulation
	if (executionMode == _PARAM_MIGRATION_ && _MIGRATION_DETAILED_OUTPUT_) {
		std::string onePop = ss.str();
		
		assert(subPopulations.size() > 0);
		
		for (std::size_t i = 0; i < subPopulations.size() - 1; ++i) {
			ss << _MIGRATION_OUTPUT_SEPARATOR_ << onePop;
		}
	}
	
	return ss.str();
}

std::size_t Simulation::getPrecision() const {
	return precision;
}

void Simulation::update(int t) {	
	switch (executionMode) {
		case _PARAM_MUTATIONS_:
			RandomDist::multinomial(allelesCount);	
			mutatePopulation();
			break;

        case _PARAM_SELECTION_:
        	updateWithSelection();
        	break;
        	
        case _PARAM_MIGRATION_:
			updateWithMigration();
			break;
			
		case _PARAM_BOTTLENECK_:
			bottleneck(t);
			RandomDist::multinomial(allelesCount, populationSize);
			break;
        	
        case _PARAM_NONE_:
		default:
			RandomDist::multinomial(allelesCount);
			break;
	}
}

void Simulation::mutatePopulation() {
	assert(!mutationFqs.empty());
		
	// mutations
	size_t nbMarkers = alleles.front().size();
	for (size_t markerIdx = 0; markerIdx < nbMarkers; ++markerIdx) {
		
		size_t nbAlleles = allelesCount.size();
		for (size_t alleleIdx = 0; alleleIdx < nbAlleles; ++alleleIdx) {
		
			// generate certain number of mutations
			int nbMut = RandomDist::binomial(allelesCount[alleleIdx], mutationFqs[markerIdx]);
			
			// iterate over mutations
			for (int mutation = 0; mutation < nbMut; ++mutation) {
										
				// generate the target mutation
				double mut = RandomDist::uniformDoubleSingle(0.0, 1.0);
				Nucl::Nucleotide target = Nucl::Nucleotide::N;
				
				double pCount = 0.0;
				for (int i = 0; i < Nucl::Nucleotide::N; ++i) {
					pCount += mutationTable[Nucl::fromChar.at(alleles[alleleIdx][markerIdx])][i];
					if (mut <= pCount) {
						target = (Nucl::Nucleotide) i;
						break;
					}							
				}
				
				if (target == Nucl::Nucleotide::N) {
					std::cerr << _ERROR_MUTATION_TARGET_UNFINDABLE_ << std::endl;
					throw _ERROR_MUTATION_TARGET_UNFINDABLE_;
				}
				
				// create new mutated allele
				std::string newAllele = std::string(alleles[alleleIdx]);
				newAllele[markerIdx] = Nucl::toChar[target];
				
				// remove original allele
				allelesCount[alleleIdx]--;
				
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
				std::cout << "Mutated a " << alleles[alleleIdx] << " to a " << newAllele << std::endl;
			}
		}
	}
}

void Simulation::bottleneck(int simulationTime) {
	if (simulationTime == bottleneckStart) {
		populationSize /= popReduction;
	} else if (simulationTime == bottleneckEnd) {
		populationSize *= popReduction;
	}
}

void Simulation::updateWithSelection() {
	// genetic drift
	int nParent = populationSize;
	int nOffspring = 0;
	double nParentCorrection = 0.0;

	assert(alleles.size() == allelesCount.size());
	assert(alleles.size() == selectionFqs.size());

	// calculation of the "corrective factor" needed to adjust the 
	// allele's frequency with the selection factor
	for (size_t i(0); i < allelesCount.size(); ++i) {
		nParentCorrection += allelesCount[i] * selectionFqs[i];
	}
	
	for (size_t i(0); i < allelesCount.size(); ++i) {
		auto& count = allelesCount[i];
		
		// remaining parent population size should be 0 or more	
		assert(nParent >= 0);
		if (nParent == 0) {
			// the only way for the parent population to be 0 is if the 
			// current allele is not present anymore (wiped out)
			assert(count == 0);
			
			// if the parent population is empty, there are no more 
			// alleles to be distributed, thus we can exit the loop
			break;
		} 
		
		// generate new allele copy number including selection frequency
		double adjustedPopulation = nParent + nParentCorrection;
		assert(adjustedPopulation != 0.0);
		
		double p = count * (1 + selectionFqs[i]) / adjustedPopulation;
		
		// reduce residual "gene pool"
		nParent -= count;
		nParentCorrection -= count * selectionFqs[i];
		
		// generate new number of allele copies in population
        count = RandomDist::binomial(populationSize - nOffspring, p);
		
		// increase offspring population size
		nOffspring += count;
	}
	
	// make sure the loop ran successfully
	assert(nParent == 0);
	assert(nOffspring == populationSize);
}

void Simulation::updateWithMigration() {
	// big container for all movements
	std::vector< std::vector< std::vector<unsigned int> > > exchange;
	
	for (int i = 0; i < (int) subPopulations.size(); ++i) {
		// exchange for current subpopulation
		std::vector< std::vector<unsigned int> > subExchange(subPopulations.size(), std::vector<unsigned int>(subPopulations[i].size(), 0));
		
		// new values for population that will go
		int gone = 0;
		for (int j = 0; j < (int) subExchange.size(); ++j) {
			gone += migrationRates[i][j];
			subExchange[j] = RandomDist::multinomialByValue(subPopulations[i], migrationRates[i][j]);
		}
	
		// new values for population that stays
		subPopulations[i] = RandomDist::multinomialByValue(subPopulations[i], subPopulationSizes[i] - gone);

		exchange.push_back(subExchange);
	}
	
	int nAlleles = subPopulations.front().size();
	
	// assign exchanges to target populations
	for (int i = 0; i < (int) exchange.size(); ++i) {
		for (int j = 0; j < (int) exchange[i].size(); ++j) {
			for (int k = 0; k < nAlleles; ++k) {
				subPopulations[j][k] += exchange[i][j][k];
			}
		}
	}	
}

int Simulation::getPopSize() const {
	return populationSize;
}

unsigned int Simulation:: subPopulationSize (std::vector<unsigned int> sub ) {

	unsigned int size(0);

	for (auto pop : sub) {

		size+=pop;

	}
	return size ;
}


std::vector< std::vector<unsigned int> > Simulation:: getSubPop() {

	return subPopulations;
}


std::vector<unsigned int> Simulation:: getAlleleCount(){

	return allelesCount;
}
