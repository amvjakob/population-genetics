#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>
#include <list>

enum Input {
	Generation, Replicas, Sites, MigrationA, MigrationT, MigrationC, MigrationG,NoInput
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
	 * Initialises populationSize, numberGenerations, replicates, numberAlleles to 0
	 * Initialises the name of the input files

	 * */
	Data ();

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
	void collectDataFile(std::ifstream& file);

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
	 * \return populationSize, an int
	 * */
	int getPopSize() const;

	/** \brief Getter of the number of generations
	 *
	 * \return numberGenerations, an int
	 * */
	double getGenerations() const;

	/** \brief Getter of the number of alleles
	 *
	 * \return numberAlleles, an int
	 * */
	int getNumberAlleles() const;

	/** \brief Getter of the number of replicates
	 *
	 * \return replicates, an int
	 * */
	double getReplicates() const;

	/** \brief Getter of the marker sites
	 *
	 * \return markerSites, a list of doubles
	 * */
	std::list<double> getMarkerSites() const;

	/** \brief Collects information related to the alleles
	 *
	 * Calculates the number of alleles and their frequencies
	 *
	 * */
	void countAlleles();

	/** \brief Getter of the nucleotide migration probabilities
	 *
	 * \return migrations, a list of double
	 * */
	std::list <double> getMigrations() const;

	/** \brief Utility function to read data from the user input file
   *
	 * */
	double comparison(double,std::string,std::string);

  std::list<double> compar_son(double,std::string,std::string);

  Input resolveInput (std::string );

  void setMigrations(std::list<double>);

	private:

	//!< Name of the user input file, a string
	std::string dataName;

	//!< Name of the fasta file, a string
	std::string fastaName;

	//!< Size of the population, an int
	int populationSize;

	//!< Number of generations (number of simulation steps), an int
	double numberGenerations;

	//!< Vector of double containing the allele frequencies of the fasta file
	std::vector<double> alleleFq;

	//!< List of double containing the user marker sites
	std::list<double> markerSites;

	//!< Number of replicates of the simulation, an int
	double replicates;

	//!< List of strings containing the allele sequences of all the individuals of the simulation
	std::list <std::string> sequences;

	//!< Number of alleles, an int
	int numberAlleles;

	//!< List of double containing the migration probabilities of respectively nucleotides A, T, C, G
	std::list <double> migrations;

};

#endif
