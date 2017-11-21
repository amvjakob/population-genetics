#ifndef ALLELE_H
#define ALLELE_H

#include <map>
#include <vector>
#include "Random.hpp"

typedef enum Nucleotide { A, C, G, T, N } Nucleotide;

class Allele {
public:

	Allele(const std::string&); //!< Constructor
	Allele(const std::string&, const std::string&); //!< Constructor
	~Allele() = default; //!< Destructor

	/**
	* @return int Allele
	*/
	std::string getIdentifier() const;

	/**
	*@param int const id
	*/
	void setIdentifier(const std::string&);
	
	/**
	 * 
	 * \return The size of the identifier
	 * */
	std::size_t size() const;
	
	const std::vector<Nucleotide>& getSequence() const;
	
	
	bool operator<(const Allele&) const;
	
	std::string getRandomGenotype(const std::vector<int>&, RandomDist&) const;
	
	static const std::map<char, Nucleotide> charToNucl;
	static const char nuclToChar[6];

private:
	std::string identifier;//!< Identifier specific to the allele
	
	std::vector<Nucleotide> sequence;
};

#endif
