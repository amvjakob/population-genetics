#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <map>
#include <array>
#include <unordered_map>
#include "Allele.hpp"
#include "Globals.hpp"

/** \brief Class representing a Simulation
 *
 * In a simulation, a population of N individuals evolves during T time
 * steps, reproducing and thus sharing a combination of alleles with the
 * future generation. We aspire to simulate that.
 *
 * */
class Simulation {

public:

	Simulation() = default;
	
	Simulation(const Simulation& other);
	
	Simulation& operator=(const Simulation& other);

	/** \brief Simulation constructor
	 *
	 * Initialises a new population genetics simulation.
	 *
	 * \param allelesCount			A vector of initial allele counts
	 * */
	Simulation(std::vector<unsigned int> allelesCount);

	/** \brief Simulation constructor
	 *
	 * Initialises a new population genetics simulation.
	 *
	 * \param alleles				List of alleles in the population
	 * \param allelesCount			Number of each allele in the population (common index with \p alleles)
	 * \param mutationFqs			Marker-specific mutation rates, in order
	 * \param nuclMutationProbs		Array of nucleotide mutation probabilities, according to one of the 3 models
	 * */			
	Simulation(const std::vector<std::string>& alleles,
				const std::vector<unsigned int>& allelesCount,
				const std::vector<double>& mutationFqs,
				const std::array< std::array<double, Nucleotide::N>, Nucleotide::N >& nuclMutationProbs);
	
	/** \brief Simulation constructor
	 *
	 * Initialises a new population genetics simulation.
	 *
	 * \param alleles				List of alleles in the population
	 * \param allelesCount			Number of each allele in the population (common index with \p alleles)
	 * \param selectionRates		List of selection rates 
	 * */			
	Simulation(const std::vector<std::string>& alleles,
				const std::vector<unsigned int>& allelesCount,
				const std::vector<double>& selectionRates);
	
	/** \brief Simulation constructor
	 *
	 * Initialises a new population genetics simulation.
	 *
	 * \param alleles				List of alleles in the population
	 * \param subPopulations		Number of each allele in the subpopulations (common index with \p alleles)
	 * \param migrationRates		Matrix of migration rates 
	 * */		
	Simulation(const std::vector<std::string>& alleles,
				const std::vector< std::vector<unsigned int> >& subPopulations,
				const std::vector< std::vector<unsigned int> >& migrationRates);
				
	/** \brief Simulation constructor
	 *
	 * Initialises a new population genetics simulation.
	 *
	 * \param alleles				List of alleles in the population
	 * \param allelesCount			Number of each allele in the population (common index with \p alleles)
	 * */
	Simulation(const std::vector<std::string>& alleles,
				const std::vector<unsigned int>& allelesCount);
				
				

	/** \brief Get the alleles in the population
	 *
	 * Get the allele distribution in the population at any given simulation step.
	 * The values represent the id's of the alleles in the population.
	 *
	 * \return A constant reference on the alleles in the population
	 *
	 * */
	const std::vector<std::string>& getAlleles() const;
	
	/** \brief Get the number of alleles in the population
	 *
	 * Get the allele distribution in the population at any given simulation step.
	 * The values represent the number of alleles in the population.
	 *
	 * \return A constant reference on the list of number alleles in the population
	 *
	 * */
	const std::vector<unsigned int>& getAllelesCount() const;

	/** \brief Utility function to format the allele numbers to frequencies for the output
	 *
	 * \return A string containing the allele frequencies at the current
	 * time step. with the following format: 0.50|...|0.25
	 * */
	std::string getAlleleFqsForOutput() const;

	/** \brief Utility function to get and format the allele identifiers
	 *
	 * \return A string containing the allele identifiers with the following
	 * format: id1|id2|..|idn-1|idn|
	 * */
    
    
    std::string getMigAlleleFqsForOutput() ;

    
	std::string getAlleleStrings() const;


	/** \brief Update the Simulation by one step
	 *
	 * "Creates" a new population of N individuals, choosing the alleles
	 * from the parent generation using a multinomial distribution.
	 *
	 * */
	void update(int t);
	
	/** \brief Get the output precision for the frequencies 
	 * */
	std::size_t getPrecision() const;
	
	/** \brief Bottleneck effet
	 *
	 * Creates a time-dependent population size
	 * 
	 * \param the time of the simulation, an int
	 * */
	void bottleneck(int simulationTime);

    /** \brief Update the Simulation by one step
	 *
	 * "Creates" a new population of N individuals, choosing the alleles
	 * from the parent generation using a multinomial distribution.
	 * Add the selection frequency to each allele.
	 *
	 * */
	void updateWithSelection();
	
	/** \brief Update the Simulation by one step
	 *
	 * */
	void updateWithMigration();
	
	/** \brief Update the Simulation by one step
	 * 
	 * */
	void updateWithMutations();

private:

	/** \brief Calculate ouput constants
	 * 
	 * Caculates the needed precision and additional spaces for the output to be aligned
	 * 
	 * */
	void calcOutputConstants();
	
	
	/** \brief Generates mutations in the current population
	 * 
	 * Mutates nucleotides from the marker sequences depending on model
	 * 
	 * */
	void mutatePopulation();
	
	
	//!< Execution mode
	int executionMode;
	
	
	//!< Size of the population
	int populationSize;

	//!< List of alleles of the current simulation
	std::vector<std::string> alleles;
	
	//!< Count of the alleles in the current simulation
	std::vector<unsigned int> allelesCount;
	
	
	//!< List of marker-specifix mutation frequencies
	std::vector<double> mutationFqs;
	
	//!< Mutation rates for every nucleotide to every nucleotide
	std::array< std::array<double, N>, N > mutationTable;
	
	
	//!< List of selections frequencies of each alleles
	std::vector<double> selectionFqs;
	
	
	//!< Table containing the sub-populations
    std::vector< std::vector<unsigned int> > subPopulations;
    
    //!< Table of subpopulation sizes
    std::vector<std::size_t> subPopulationSizes;
    
    //!< Table containing migration rates for each sub group
    std::vector< std::vector<unsigned int> > migrationRates;
	
	
	//!< Precision for output
	std::size_t precision;
	
	//!< Additional spaces for correct output format
	std::size_t additionalSpaces;
};


#endif
