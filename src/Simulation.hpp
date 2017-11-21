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
	 * \param alleleFq				A vector of initial allele frequencies
	 * */
	Simulation(int populationSize, std::vector<double> alleleFq);

	/** \brief Simulation constructor
	 *
	 * Initialises a new population genetics simulation.
	 *
	 * \param alleles				A map of alleles vs number in the population
	 * */
	Simulation(const std::unordered_map<std::string, int>& alleles);

	/** \brief Get the allele distribution in the population
	 *
	 * Get the allele distribution in the population at any given simulation step.
	 * The returned map's indices indicate the index of the allele,
	 * the map's values represent the number of this allele in the population.
	 *
	 * \return A constant reference on the alleles in the population
	 *
	 * */
	const std::unordered_map<std::string, int>& getAlleles() const;

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


private:

	/** \brief Calculate ouput constants
	 * 
	 * Caculates the needed precision and additional spaces for the output to be aligned
	 * 
	 * */
	void calcOutputConstants();
	
	
	static std::array< std::array<double, N >, N > mutationTable;


	//!< Size of the population
	int populationSize;

	//!< List of alleles of the current simulation
	std::unordered_map<std::string, int> alleles;
	
	//!< Execution mode
	const int executionMode = _PARAM_MUTATIONS_;
	
	//!< Precision for output
	std::size_t precision;
	
	//!< Additional spaces for correct output format
	std::size_t additionalSpaces;
};


#endif
