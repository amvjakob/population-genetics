#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>
#include <list>

/** Input enum 
 *  Handles the different types of data that are read 
 */
enum Input {
	Generation, Replicas, Sites, Mutations, NoInput
};


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
	double getGenerations() const;

	/** \brief Getter of the number of alleles
	 *
	 * 	\return numberAlleles, an int
	 * */
	int getNumberAlleles() const;

	/** \brief Getter of the number of replicates
	 *
	 * 	\return replicates, an int
	 * */
	double getReplicates() const;

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
	 * 	\param the size of the string to read, an int
	 * 	\param the string to read from
	 * 	\param the string of the data wanted/to be read
	 * 
	 * 	\return the data read, a double
	 * */
	double extractDouble(std::string) const;
	
	/** \brief Utility function to read data from the user input file
	 *
	 * 	\param the size of the string to read, an int
	 * 	\param the string to read from
	 * 	\param the string of the data wanted/to be read
	 * 
	 * 	\return the data read, a vector of double
	 * */
    std::vector<double> extractVec(std::string) const;
	
	/** \brief Utility function to convert a string into the input enum type
	 *
	 * 	\param the string to be converted
	 * 
	 * 	\return the enum value corresponding to a given string, an Input
	 * */
    Input resolveInput(std::string);
	
	/** \brief Setter of the vector of mutations
	 * 
	 * 	\param a vector of double, the new value of mutations
	 * */
    void setMutations(std::vector<double>&);
	
	/** \brief Getter of the vector of the allele frequencies
	 * 
	 * 	\return a vector of double, the allele frequencies
	 * */
	const std::vector<double>& getAllelesFq() const;
	
	
private:

	//!< Name of the user input file, a string
	std::string inputName;

	//!< Name of the fasta file, a string
	std::string fastaName;

	//!< Size of the population, an int
	int populationSize;

	//!< Number of generations (number of simulation steps), an int
	double numberGenerations;
	
	//!< Number of replicates of the simulation, an int
	double replicates;
	
	//!< Number of alleles, an int
	int numberAlleles;
	
	//!< Vector of double containing the allele frequencies of the fasta file
	std::vector<double> alleleFq;

	//!< Vector of double containing the user marker sites
	std::vector<double> markerSites;

	//!< List of strings containing the allele sequences of all the individuals of the simulation
	std::list<std::string> sequences;

	//!< Vector of double containing the mutations probabilities of respectively nucleotides A, T, C, G
	std::vector<double> mutations;

};

#endif
