#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <algorithm>
#include "Data.hpp"
#include "Random.hpp"

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
	bottleneckStart = 0;
	bottleneckEnd = 0;
	popReduction = 1;
	
	executionMode = _PARAM_NONE_;
	
	// mutation
	mutationModel = _MUTATION_MODEL_NONE_;

	//migration model
	migrationModel = _MIGRATION_MODEL_NONE_;

	kimuraDelta = 0.0;
	
	assert(allelesCount.empty());
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
		cerr << _ERROR_INPUT_UNREADABLE_ << endl;
		throw _ERROR_INPUT_UNREADABLE_;
	}

	if (fastaFile.is_open()) {
		collectFastaFile(fastaFile);
	} else {
		cerr << _ERROR_FASTA_UNREADABLE_ << endl;
		throw _ERROR_FASTA_UNREADABLE_;
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
		if (line[0] == _INPUT_COMMENT_) continue;
		
		stringstream ss(line);
		getline(ss, key, _INPUT_DECLARATION_);

        switch (str2int(key.c_str())) {
			case str2int(_INPUT_KEY_GENERATIONS_):
                extractValue<int>(numberGenerations, line, strToInt);
                break;
                
			case str2int(_INPUT_KEY_REPLICAS_):
				extractValue<int>(replicates, line, strToInt);
				break;

            case str2int(_INPUT_KEY_MARKER_SITES_):
				extractValues<int>(markerSites, line, strToInt);
				break;
				
			case str2int(_INPUT_KEY_MODE_):
				extractValue<int>(executionMode, line, strToInt);
				break;

			case str2int(_INPUT_KEY_MIGRATION_MODEL_):
				extractValue<int>(migrationModel, line, strToInt);
				break;

			case str2int(_INPUT_KEY_MIGRATION_MODE_):
				extractValue<int>(migrationMode, line, strToInt);
				break;

            case str2int(_INPUT_KEY_MIGRATION_RATES_):
                extractValues<int>(migrationRates, line, strToInt);
                break;

            case str2int(_INPUT_KEY_MUTATION_RATES_):
				extractValues<double>(mutations, line, strToDouble);
                break;
                
            case str2int(_INPUT_KEY_MUTATION_KIMURA_):
				extractValue<double>(kimuraDelta, line, strToDouble);
				break;
				
			case str2int(_INPUT_KEY_MUTATION_FELSENSTEIN_):
				extractValues<double>(felsensteinConstants, line, strToDouble);
				break;

			case str2int(_INPUT_KEY_SELECTION_RATES_):
				extractValues<double>(selections, line, strToDouble);
				break;
				
			case str2int(_INPUT_KEY_BOTTLENECK_POPULATION_REDUCTION_):
				extractValue<double>(popReduction, line, strToDouble);
				break;
				
			case str2int(_INPUT_KEY_BOTTLENECK_START_TIME_):
				extractValue<int>(bottleneckStart, line, strToInt);
				break;
				
			case str2int(_INPUT_KEY_BOTTLENECK_END_TIME_):
				extractValue<int>(bottleneckEnd, line, strToInt);
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
		} else if (!felsensteinConstants.empty() && felsensteinConstants.size() == (int) Nucl::Nucleotide::N) {
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
				for (auto& c : felsensteinConstants) c += (1.0 - sum) / (Nucl::Nucleotide::N);
			}
			
			// set mutation model
			if (!(sum > 1.0)) {
				mutationModel = _MUTATION_MODEL_FELSENSTEIN_;
			}
		} 
	} 
	
	// additional checks needed!
	// check for data completeness - right params, enough params
}



void Data::collectFastaFile(ifstream& file) {
	string line;

	while (getline(file, line)) {
		
		if (line[0] == _FASTA_COMMENT_) {
			//increasing population thanks to the header
			++populationSize;
			continue;
		}

		string seq;

		for (auto& marker : markerSites) {
			char c = line[marker - 1];
			
			if (std::string("ACTG").find(c) != string::npos) {
				seq += line[marker - 1];
			}
			else {
				// if we have an unknown nucleotide, generate a valid one randomly
				seq += Nucl::toChar[RandomDist::uniformIntSingle(0, 3)];
			}
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

int Data::getReplicates() const {
	return replicates;
}

void Data::countAlleles() {
	list<string> sequencesSorted = sequences;

	//list container functions that allow to sort and remove the duplicates
	sequencesSorted.sort();
	sequencesSorted.unique();

	for (auto& seqS : sequencesSorted) {
		auto count = (unsigned int) count_if(sequences.begin(), sequences.end(), [&](string allele) {
				return allele == seqS;
		});

		allelesCount.push_back(count);
	}
}

const vector<double>& Data::getMutations() const {
	return mutations;
}

const std::vector<double>& Data::getSelections() const {
	return selections;
}

const std::vector<unsigned int>& Data::getAllelesCount() const {
	return allelesCount;
}

const std::list<std::string>& Data::getSequences() const {
	return sequences;
}

std::vector<std::string> Data::getUniqueSequences() const {
	std::list<std::string> uniqueSequences = sequences;
	
	uniqueSequences.sort();
	uniqueSequences.unique();
	
	std::vector<std::string> res;
	
	for (auto& seq : uniqueSequences)
		res.push_back(seq);
	
	return res;
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

int Data::getMigrationModel() const {
	return migrationModel;
}

int Data::getMigrationMode() const {
	return migrationMode;
}

const std::vector<int>& Data::getMigrations() const {
    return migrationRates;
}

double Data::getPopReduction() const {
	return popReduction;
}

int Data::getBottleneckStart() const {
	return bottleneckStart;
}

int Data::getBottleneckEnd() const {
	return bottleneckEnd;
}

const std::vector<int>& Data::getMarkerSites() const {
	return markerSites;
}
