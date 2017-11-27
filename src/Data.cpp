#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include "Data.hpp"
#include "Allele.hpp"
#include "Globals.hpp"

using namespace std;

Data::Data(string input, string fasta) {
	construct(input, fasta);	
}

Data::Data() {
	string input;
	string fasta;
	
	
	cout << "Please enter the path of your input file: " << endl;
	cin >> input;

	cout << "Please enter the path of your fasta file: " << endl;
	cin >> fasta;
	
	construct(input, fasta);
}

void Data::construct(string input, string fasta) {
	inputName = input;
	fastaName = fasta;
	
	populationSize = 0;
	numberGenerations = 0;
	replicates = 0;
	numberAlleles = 0;
	
	executionMode = _PARAM_NONE_;
	
	// mutation
	mutationModel = _MUTATION_MODEL_NONE_;
	kimuraDelta = 0.0;
	
	assert(alleleFq.empty());
	assert(markerSites.empty());
	assert(sequences.empty());
	assert(mutations.empty());
	assert(selections.empty());
}

void Data::collectAll() {	
	ifstream dataFile;
	dataFile.open(inputName, ifstream::in);
	
	ifstream fastaFile;
	fastaFile.open(fastaName, ifstream::in);

	if (dataFile.is_open()) {
		collectUserFile(dataFile);
	} else {
		cerr << "Ouverture de l'input impossible" << endl;
		throw 10;
	}

	if (fastaFile.is_open()) {
		collectFastaFile(fastaFile);
	} else {
		cerr << "Ouverture  du fasta impossible" << endl;
		cerr << fastaName << endl;
		throw 11;
	}
}

void Data::collectUserFile(ifstream& file) {
	string line, key;

	while (getline(file, line)) {
		
		// remove whitespace
		line.erase(
			remove_if(line.begin(), line.end(), ::isspace), 
			line.end()
		);
		
		// lines starting with # are comments
		if (line[0] == '#') continue;
		
		stringstream ss(line);
		getline(ss, key, '=');

        switch (str2int(key.c_str())) {
			case str2int(_INPUT_KEY_GENERATIONS_):
                numberGenerations = extractInt(line);
                break;
                
			case str2int(_INPUT_KEY_REPLICAS_):
				replicates = extractInt(line);
				break;

            case str2int(_INPUT_KEY_MARKER_SITES_):
				markerSites = extractVec(line);
				break;
				
			case str2int(_INPUT_KEY_MODE_):
				executionMode = extractInt(line);
				break;

			case str2int(_INPUT_KEY_MUTATION_RATES_):
				mutations = extractVec(line);
                break;
                
            case str2int(_INPUT_KEY_MUTATION_KIMURA_):
				kimuraDelta = extractDouble(line);
				break;
				
			case str2int(_INPUT_KEY_MUTATION_FELSENSTEIN_):
				felsensteinConstants = extractVec(line);
				break;

			case str2int(_INPUT_KEY_SELECTION_RATES_):
				selections = extractVec(line);
				break;

            default:
				break;
        }
	}
	
	// the whole file has been read
	
	// set mutation model, if mutation_mode
	if (executionMode == _PARAM_MUTATIONS_) {
		// default is cantor
		mutationModel = _MUTATION_MODEL_CANTOR_;
		
		if (kimuraDelta >= 1.0/3.0 && kimuraDelta <= 1.0) {
			mutationModel = _MUTATION_MODEL_KIMURA_;
		} else if (!felsensteinConstants.empty() && felsensteinConstants.size() == (int) Nucleotide::N) {
			// check for correct constants
			double sum = 0.0;
			for (auto& c : felsensteinConstants) {				
				// c can not be negative
				c = std::abs(c);
				
				// c can not be equal to 1
				sum += c;
			}
			
			// adjust terms
			if (sum < 1.0) {
				for (auto& c : felsensteinConstants) c += (1.0 - sum) / (Nucleotide::N);
			}
			
			// set mutation model
			if (!(sum > 1.0)) {
				mutationModel = _MUTATION_MODEL_FELSENSTEIN_;
			}
		} 
	}
}



void Data::collectFastaFile(ifstream& file) {
	string line;

	while (getline(file, line)) {
		
		if (line[0] == '>') {
			++populationSize;
			continue;
		}

		string seq;

		for (auto& marker : markerSites) {
			seq += line[marker - 1];
		}

		sequences.push_back(seq);
	}

	countAlleles();
}

int Data::getPopSize() const {
	return populationSize;
}

int Data::getGenerations() const {
	return numberGenerations;
}

int Data::getNumberAlleles() const {
	return numberAlleles;
}

int Data::getReplicates() const {
	return replicates;
}

const std::vector<double>& Data::getMarkerSites() const {
	return markerSites;
}

void Data::countAlleles() {
	list<string> sequencesSorted = sequences;
	
	sequencesSorted.sort();
	sequencesSorted.unique();
	
	numberAlleles = sequencesSorted.size();

	for (auto& seqS : sequencesSorted) {
		int count = count_if(sequences.begin(), sequences.end(), [&](string allele) {
				return allele == seqS;
		});

		alleleFq.push_back(count * 1.0 / populationSize);
	}
}

/*void Data::setSelectionsToAlleles() {
	
	for(size_t i(0); i < selections.size(); ++i) {
		double allWithSelec = allelesNum[i]*selections[i];
		total += allWithSelec;
	}

	double ajustedPopulation = populationSize + total;

	for (size_t i(0); i < selections.size(); ++i) {
		double selectionAllele = allelesNum[i]*selections[i];

		alleleFq.push_back(selectionAllele / ajustedPopulation);
	}
}*/

const vector<double>& Data::getMutations() const {
	return mutations;
}

const std::vector<double>& Data::getSelections() const {
	return selections;
}

int Data::extractInt(string line) const {
	string key, strValue;
	
	int value = 0;
    
    stringstream ss(line);
    getline(ss, key, '=');
	getline(ss, strValue);
	
	try {
		value = stoi(strValue);
	} catch (std::invalid_argument& e) {	
		cerr << e.what() << endl;
	}

	return value;
}

double Data::extractDouble(string line) const {
	string key, strValue;
	
	double value = 0;
    
    stringstream ss(line);
    getline(ss, key, '=');
	getline(ss, strValue);
	
	try {
		value = stod(strValue);
	} catch (std::invalid_argument& e) {	
		cerr << e.what() << endl;
	}

	return value;
}

std::vector<double> Data::extractVec(string line) const {
    string key, strValue;
    
    vector<double> values;

	stringstream ss(line);
	getline(ss, key, '=');

	while (getline(ss, strValue, '|')) {
		try {
			values.push_back(stod(strValue));
		} catch (std::invalid_argument& e) {
			cerr << e.what() << endl;
		}
	}

    return values;
}

void Data::setMutations(std::vector<double>& list) {
    for (auto mig : list) {
        mutations.push_back(mig);
    }
}

const std::vector<double>& Data::getAlleleFqs() const {
	return alleleFq;
}

const std::list<std::string>& Data::getSequences() const {
	return sequences;
}

int Data::getExecutionMode() const {
	return executionMode;
}

int Data::getMutationModel() const {
	return mutationModel;
}

double Data::getKimuraDelta() const {
	return kimuraDelta;
}

const std::vector<double>& Data::getFelsensteinConstants() const {
	return felsensteinConstants;
}
