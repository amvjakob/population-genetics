#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>
#include <list>


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
	 * Initialises populationSize, numberGenerations, replicates, numberAlleles to 0
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

	/** \brief Getter of the marker sites
	 *
	 * 	\return markerSites, a vector of doubles
	 * */
	const std::vector<double>& getMarkerSites() const;

	/** \brief Collects information related to the alleles
	 *
	 * Calculates the number of alleles and their frequencies
	 *
	 * */
	void countAlleles();

	/** \brief Getter of the sites mutations probabilities
	 *
	 * 	\return migrations, a vector of double
	 * */
	const std::vector<double>& getMutations() const;

	/** \brief Utility function to read data from the user input file
	 *
	 * 	\param the string to read from
	 * 
	 * 	\return the data read, an int
	 * */
	int extractInt(std::string) const;
	
	/** \brief Utility function to read data from the user input file
	 *
	 * 	\param the string to read from
	 * 
	 * 	\return the data read, a double
	 * */
	double extractDouble(std::string) const;
	
	/** \brief Utility function to read data from the user input file
	 *
	 * 	\param the string to read from
	 * 
	 * 	\return the data read, a vector of double
	 * */
    std::vector<double> extractVec(std::string) const;
	
	/** \brief Setter of the vector of mutations
	 * 
	 * 	\param a vector of double, the new value of mutations
	 * */
    void setMutations(std::vector<double>&);
	
	/** \brief Getter of the vector of the allele frequencies
	 * 
	 * 	\return a vector of double, the allele frequencies
	 * */
	const std::vector<double>& getAlleleFqs() const;
	
	/** \brief Getter of the list of allele sequences
	 * 
	 * 	\return a list of strings, the allele sequences
	 * */
	const std::list<std::string>& getSequences() const;
	
	/** \brief Utility function to transfrom strings to ints for use in switch statements
	 * 
	 * */
	static constexpr unsigned int str2int(const char* str, int h = 0) {
		return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
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
	
	//!< Number of alleles, an int
	int numberAlleles;
	
	//!< Vector of double containing the allele frequencies of the fasta file
	std::vector<double> alleleFq;

	//!< Vector of double containing the user marker sites
	std::vector<double> markerSites;

	//!< List of strings containing the allele sequences of all the individuals of the simulation
	std::list<std::string> sequences;

	//!< Vector of double containing the mutations probabilities of the marker sites
	std::vector<double> mutations;
	
	//!< Execution mode (param to use)
	int executionMode;
	
	//!< Mutation model (simple, kimura, felsenstein)
	int mutationModel;
	
	//!> Kimura model
	double kimuraDelta;
	
	//!< Felsenstein model
	std::vector<double> felsensteinConstants;

};

#endif
