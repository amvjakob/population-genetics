#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include "Data.hpp"

using namespace std;

Data::Data(string input, string fasta)
  : inputName(input), fastaName(fasta), populationSize(0), 
  numberGenerations(0), replicates(0), numberAlleles(0)
{
	assert(alleleFq.empty());
	assert(markerSites.empty());
	assert(sequences.empty());
	assert(mutations.empty());
}

Data::Data() {
	string input;
	string fasta;
	
	cout << "Please enter the path of your input file: " << endl;
	cin >> input;

	cout << "Please enter the path of your fasta file: " << endl;
	cin >> fasta;
	
	Data(input, fasta);
}

void Data::collectAll() {
	ifstream dataFile(inputName, ios::in);
	ifstream fastaFile(fastaName, ios::in);

	if (dataFile.is_open()) {
		collectUserFile(dataFile);
	} else {
		cerr << "Ouverture de l'input impossible" << endl;
	}

	if (fastaFile.is_open()) {
		collectFastaFile(fastaFile);
	} else {
		cerr << "Ouverture  du fasta impossible" << endl;
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

        switch (resolveInput(key)) {
			case Generation:
                numberGenerations = extractDouble(line);
                break;
                
			case Replicas:
				replicates = extractDouble(line);
				break;

            case Sites:
				markerSites = extractVec(line);
				break;

			case Mutations:
				mutations = extractVec(line);
                break;

            case NoInput:
            default:
				break;
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

double Data::getGenerations() const {
	return numberGenerations;
}

int Data::getNumberAlleles() const {
	return numberAlleles;
}

double Data::getReplicates() const {
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

	double count(0);
	
	
	for (auto& seqS : sequencesSorted) {
		count = count_if(sequences.begin(), sequences.end(), [&](string allele) {
				return allele == seqS;
		});

		alleleFq.push_back(count / populationSize);
	}
}

const vector<double>& Data::getMutations() const {
	return mutations;
}

double Data::extractDouble(string line) const {
	string key, strValue;
	
	double value(0);
    
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

Input Data::resolveInput(std::string input) {
    if( input == "GEN" ) return Generation;
    if( input == "REP" ) return Replicas;
    if( input == "SITES")return Sites;
    if( input == "MUT" ) return Mutations;

    return NoInput ;
}

void Data::setMutations(std::vector<double>& list) {
    for (auto mig : list) {
        mutations.push_back(mig);
    }
}

const std::vector<double>& Data::getAllelesFq() const {
	return alleleFq;
}
