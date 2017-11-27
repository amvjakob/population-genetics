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

	/** \brief Simulation constructor
	 *
	 * Initialises a new population genetics simulation.
	 *
	 * \param populationSize		N, the number of individuals in the population
	 * \param allelesCount			A vector of initial allele counts
	 * */
	Simulation(int populationSize, std::vector<int> allelesCount);

	/** \brief Simulation constructor
	 *
	 * Initialises a new population genetics simulation.
	 *
	 * \param alleles				Map of alleles in the population, the mapped value represents the number of said allele in the population
	 * \param executionMode			The param to use for this simulation (mutation, migration, ...)
	 * \param mutationFqs			Marker-specific mutation rates, in order
	 * \param nuclMutationProbs		Array of nucleotide mutation probabilities, according to one of the 3 models
	 * */
	Simulation(const std::unordered_map<std::string, int>& alleles,
			const int executionMode,
			const std::vector<double>& mutationFqs, const std::array< std::array<double, Nucleotide::N>, 
			Nucleotide::N >& nuclMutationProbs, std::vector<double>& selectionRates);

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
	void update();
	
	/** \brief Get the output precision for the frequencies 
	 * */
	std::size_t getPrecision() const;
	
	/** \brief Bottleneck effet
	 *
	 * Creates a time-dependent population size
	 * 
	 * \param the time of the simulation, an int
	 * */
	void bottleneck (int simulationTime);
    
    /** \brief Assign each allele to a sub group
     *
     * Implementation of the sub group table using the alleleCount list
     *
     *Each subgroup takes a certain number of each alleles
     *
     * */
    void subPopCreation () ;
    
    /** \brief Creates migration rate to between each subgroup
     *
     * Implementation of the migration matrix
     *
     *checking that there is no negative rate
     *
     * */
    void migrationRatesCreation();

	
    /** \brief Migration effect
     *
     * Creates a migrations between subGroups
     *
     *alleles exchanges
     *
     * */
    void migrationUpdate();
    
    /** \brief Calculate subgroup size
     *
     *
     *make sure population size stays constant
     *
     * */
    double subgroupSize(std::vector< double >);
    
    /** \brief prepare migrations vectors
     *
     *intialize their size
     *
     * */
    void prepareMigrationVectors ();

    /** \brief Update the Simulation by one step
	 *
	 * "Creates" a new population of N individuals, choosing the alleles
	 * from the parent generation using a multinomial distribution.
	 * Add the selection frequency to each allele.
	 *
	 * */
	void updateWithSelection();

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
	
	
	//!< Size of the population
	int populationSize;

	//!< List of alleles of the current simulation
	std::vector<std::string> alleles;
	
	//!< Count of the alleles in the current simulation
	std::vector<unsigned int> allelesCount;
	
	//!< Execution mode
	const int executionMode;
	
	//!< List of marker-specifix mutation frequencies
	std::vector<double> mutationFqs;
	
	//!< Mutation rates for every nucleotide to every nucleotide
	std::array< std::array<double, N>, N > mutationTable;
	
	
	//!< Precision for output
	std::size_t precision;
	
	//!< Additional spaces for correct output format
	std::size_t additionalSpaces;
    
    
    
    
    
    //!< Table containing the sub populations
    std::vector<std::vector<double>> subPops;
    
    //! < Table containing migration rates for each sub group
    std::vector<std::vector<double>> migrationTable;



    	//!< list of selections frequencies of each alleles
	std::vector<double> selectionFqs;
};


#endif
