#include "Allele.hpp"

const std::map<char, Nucleotide> Allele::charToNucl = { 
	{'a', A}, {'A', A},
	{'c', C}, {'C', C},
	{'g', G}, {'G', G},
	{'t', T}, {'T', T},
	{'n', N}, {'N', N}, {'x', N}, {'*', N}
};

const char Allele::nuclToChar[6] = "ACGTN";

Allele::Allele(const std::string& id)
  : identifier(id)
{ 
	// empty sequence
}

Allele::Allele(const std::string& id, const std::string& seq)
  : identifier(id)
{ 
	// build sequence
	for (auto it = seq.begin(); it != seq.end(); ++it) {
		sequence.push_back(Allele::charToNucl.at(*it));
	}
}

std::string Allele::getIdentifier() const {
	return identifier;
}

void Allele::setIdentifier(const std::string& id) {
	identifier = id;
}

bool Allele::operator<(const Allele& other) const {
	return identifier < other.identifier;
}

std::string Allele::getRandomGenotype(const std::vector<int>&, RandomDist&) const {
	return "";
}

