#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>
#include <list>
#include <functional> 
#include <iostream>
#include <sstream>
#include "Globals.hpp"


/** \brief Class regrouping the data necessary to run a simulation
 *
 * Collecting the data requires : a user input file (.txt) and a fasta file (.fa)
 *
 * */
class Data {

public:

	/** \brief Data constructor
	 *
	 * 	\param input, the path of the input file to be read
	 * 	\param fasta, the path of the fasta file to be read
	 * 
	 * Initialises populationSize, numberGenerations, replicates, numberAlleles, bottleneckStart, bottleneckEnd to 0 and popReduction to 1
	 * Initialises the name of the data files, according to the parameters
	 * */
	Data(std::string input, std::string fasta);
	
	
	/** \brief Second Data constructor
	 *
	 * Asks the file names to the user
	 * 
	 * Call the first constructor
	 * 
	 * */
	Data();
	
	
	/** \brief Construct a new Data object
	 * */
	void construct(std::string input, std::string fasta);
	
	/** \brief Collects all the data required for the program
	 *
	 * Collects inputs from the user and the fasta files
	 * */
	void collectAll();

	/** \brief Collects data from the user file
	 *
	 * Reads the number of generations, the marker sites, the number of getReplicates
	 * Also reads the migration probabilities of each nucleotide
	 * */
	void collectUserFile(std::ifstream& file);

	/** \brief Collects data from the fasta file
	 *
	 * Calculates the number of individuals/size of the population
	 * Registers the allele sequences of all the individuals
	 * Counts the number of different alleles and their frequencies
	 * Initialises the name of the input files
	 *
	 * */
	void collectFastaFile(std::ifstream& file);

	/** \brief Getter of the size of the populationSize
	 *
	 * 	\return populationSize, an int
	 * */
	int getPopSize() const;

	/** \brief Getter of the number of generations
	 *
	 * 	\return numberGenerations, an int
	 * */
	int getGenerations() const;

	/** \brief Getter of the number of alleles
	 *
	 * 	\return numberAlleles, an int
	 * */
	int getNumberAlleles() const;

	/** \brief Getter of the number of replicates
	 *
	 * 	\return replicates, an int
	 * */
	int getReplicates() const;

	/** \brief Collects information related to the alleles
	 *
	 * Calculates the number of alleles and their frequencies
	 *
	 * */
	void countAlleles();

	/** \brief Getter of the sites mutations probabilities
	 *
	 * 	\return mutations , a vector of double
	 * */
	const std::vector<double>& getMutations() const;

	/** \brief Get selection rates
	 * */
	const std::vector<double>& getSelections() const;
	
	/** \brief Getter of the vector of the allele counts
	 * 
	 * 	\return a vector of int, the allele counts
	 * */
	const std::vector<unsigned int>& getAllelesCount() const;
	
	const std::vector<int>& getMarkerSites() const;
	
	/** \brief Getter of the list of allele sequences
	 * 
	 * 	\return a list of strings, the allele sequences
	 * */
	const std::list<std::string>& getSequences() const;
	
	/** \brief Get the different alleles in the Simulation
	 * 
	 * \return a vector of strings, the alleles
	 * */
	std::vector<std::string> getUniqueSequences() const;
	
	/** \brief Utility function to transfrom strings to ints for use in switch statements
	 * 
	 * */
	static constexpr unsigned int str2int(const char* str, int h = 0) {
		return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
	} 
	
	/** \brief Get the execution mode of a sSmulation (mutation, migration, ...)
	 * 
	 * \return An int whose meaning is defined in Globals.hpp
	 * */
	int getExecutionMode() const;
	
	/** \brief Get the mutation model to use for a Simulation (Cantor, Kimura, ...)
	 * 
	 * */
	int getMutationModel() const;
	
	/** \brief Get the value necessary to create a Kimura mutation model
	 * 
	 * \return The value of delta for a Kimura mutation model
	 * */
	double getKimuraDelta() const;
	
	/** \brief Get the values necessary to create a Felsenstein mutation model
	 * 
	 * \return The values of the constants for a Felsenstein model
	 * */
	const std::vector<double>& getFelsensteinConstants() const;
	

    /** \brief Get the migration model to use for a Simulation (comppleteGraph, Star, ...)
	 *
	 * */
    int getMigrationModel() const;

	/** \brief Get the migration mode to use for a Simulation (random, user input , ...)
	 *
	 * */
	int getMigrationMode() const;
	
	/** \brief Getter of the sites migration rates for each allele
     *
     * user supposed to know number of different allele
     *
     * each rate correspond to the amount of outgoing people
     *
     * 	\return migrations, a vector of double
     * */
    const std::vector<int>& getMigrations() const;
    
    
    /** \brief Extracts a value from a line
	 * 
	 * \param into			the variable to store the extracted value into
	 * \param line			the line to read from, a string
	 * \param extractFn		the function to call to cast the value contained in the string to the corresponding arithmetic type
	 * */
    template<typename T>
	void extractValue(T& into, std::string line, std::function< T (const std::string&) > extractFn) {
		std::string key, strValue;

		std::stringstream ss(line);
		std::getline(ss, key, _INPUT_DECLARATION_); // separate key from value
		std::getline(ss, strValue); // store value in strValue
		
		try {
			// cast value
			into = extractFn(strValue);
			
		} catch (std::invalid_argument& e) {	
			std::cerr << e.what() << std::endl;
			throw e.what();
		}
	}
	
	/** \brief Extracts a vector of values from a line
	 * 
	 * \param into			the variable to store the extracted values into
	 * \param line			the line to read from, a string
	 * \param extractFn		the function to call to cast the value contained in the string to the corresponding arithmetic type
	 * */
	template<typename T>
	void extractValues(std::vector<T>& into, std::string line, std::function< T (const std::string&) > extractFn) {
		std::string key, strValue;

		std::stringstream ss(line);
		std::getline(ss, key, _INPUT_DECLARATION_); // separate key from value

		while (std::getline(ss, strValue, _INPUT_SEPARATOR_)) {
			try {
				//add elements found between each separators
				T extractedValue = extractFn(strValue);
				into.push_back(extractedValue);

			} catch (std::invalid_argument& e) {
				std::cerr << e.what() << std::endl;
				throw e.what();
			}
		}
	}
	
	/** \brief Getter of the reduction of the population size factor during the bottleneck
	 *
	 * 	\return popReduction, a double
	 * */
	 double getPopReduction() const;
	 
	 /** \brief Getter of the start of the bottleneck
	 *
	 * 	\return bottleneckStart, an int
	 * */
	 int getBottleneckStart() const;
	 
	 /** \brief Getter of the end of the bottleneck
	 *
	 * 	\return bottleneckEnd, an int
	 * */
	 int getBottleneckEnd() const;
	 
	
private:

	//!< Name of the user input file, a string
	std::string inputName;

	//!< Name of the fasta file, a string
	std::string fastaName;

	//!< Size of the population, an int
	int populationSize;

	//!< Number of generations (number of simulation steps), an int
	int numberGenerations;
	
	//!< Number of replicates of the simulation, an int
	int replicates;
	
	//!< Vector of double containing the allele frequencies of the fasta file
	std::vector<unsigned int> allelesCount;

	//!< Vector of double containing the user marker sites
	std::vector<int> markerSites;

	//!< List of strings containing the allele sequences of all the individuals of the simulation
	std::list<std::string> sequences;


	//!< Execution mode (param to use)
	int executionMode;
	

	//!< Vector of double containing the mutations probabilities of the marker sites
	std::vector<double> mutations;
	
	//!< Mutation model (simple, kimura, felsenstein)
	int mutationModel;
	
	//!< Kimura model
	double kimuraDelta;
	
	//!< Felsenstein model
	std::vector<double> felsensteinConstants;


	//!< Migration  mode (user input , random )
	int migrationMode;

    //!< Migration  model (simple, kimura, felsenstein)
    int migrationModel;
    
    //!< Migration rates
    std::vector<int> migrationRates;


	//!< Vector of double containing the selection probabilities of the alleles
	std::vector<double> selections;

	
	//!< Bottleneck population reduction factor
	double popReduction;
	
	//!< Bottleneck start time
	int bottleneckStart;
	
	//!< Bottleneck stop time
	int bottleneckEnd;
};

#endif
