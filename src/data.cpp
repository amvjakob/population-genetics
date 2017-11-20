#include "data.hpp"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cassert>

using namespace std;

Data::Data(string input, string fasta) :  inputName (input),fastaName(fasta),populationSize(0), numberGenerations (0), replicates (0),numberAlleles(0) {

	assert(alleleFq.empty());
	assert(markerSites.empty());
	assert (sequences.empty());
	assert (mutations.empty());
	

}

Data::Data() {
	
	string input;
	string fasta;
	
	cout << "Please enter the path of your input file: " << endl;
	cin >> input;

	cout << "Please enter the path of your fasta file: " << endl;
	cin >> fasta;
	
	Data (input, fasta);

}
	

void Data::collectAll() {

	ifstream dataFile (inputName, ios::in);
	ifstream fastaFile (fastaName, ios::in);

	if (dataFile.is_open()) {

		collectUserFile(dataFile);
	} else {
		cout << "Ouverture  de l'input impossible" << endl;
	}

	if (fastaFile.is_open()) {

		collectFastaFile(fastaFile);

	} else {
		cout << "Ouverture  du fasta impossible" << endl;
	}

}

void Data::collectUserFile(ifstream& file) {


  string line;
	string key;


	while (getline(file,line)) {
		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
		if (line[0] == '#') continue;
		stringstream ss(line);
		getline(ss, key, '=');

        switch (resolveInput(key))
        {

               case Generation :
               {

                numberGenerations=comparison(3,"GEN",line);

            }

                case Replicas :
                {
                    replicates=comparison(3,"REP",line);
                }

                case Sites :
                {
                    markerSites=compar_son(5,"SITES",line);
                }


                case Mutations :
                {
                    setMutations(compar_son(3,"MUT",line));
                }

                 case NoInput :
                 {
                     break;

                 }


        }

        }

	}



void Data::collectFastaFile(ifstream& file) {

	string line;

	while (getline(file,line)) {
		if (line[0] == '>') {
			++populationSize;
			continue;
		}

		string seq;

		for ( auto marker : markerSites ) {
			seq += line [marker-1];
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

double  Data::getReplicates() const {
	return replicates;
}

std::vector<double> Data::getMarkerSites() const {
	return markerSites;
}

void Data::countAlleles() {

	list <string> sequencesSorted = sequences;
	sequencesSorted.sort();
	sequencesSorted.unique();
	numberAlleles = sequencesSorted.size();

	double count (0);

	for (auto seqS : sequencesSorted) {
		count = 0;
		for (auto seq : sequences) {
			if (seq == seqS) {
				++count;
			}
		}

		alleleFq.push_back(count/populationSize);
	}

}

vector <double> Data::getMutations() const {
	return mutations;
}

double Data:: comparison(int size,string s2,string l)
{
	string value ;
    string key;
	double v(0);
          stringstream ss(l);
         getline(ss, key, '=');
   if (key.compare(0,size,s2)==0)
   {
       std::getline(ss, value);
       try
       {
           v = stoi(value);
       }

       catch (std::invalid_argument& e)
       {}
   }

          return v;

      }




std::vector<double> Data::compar_son(int size,std::string s2,std::string l)
{
    string key;
    string value ;
    std::vector<double> v;

        stringstream ss(l);
        getline(ss, key, '=');

    if (key.compare(0,size,s2)==0) {

        if (size == 5) {

            while (getline(ss, value, '|'))
            {
              try {
                  v.push_back(stod(value));
              }

              catch (std::invalid_argument& e)
              {}
            }
        } else {

            std::getline(ss,value);
            try {
                v.push_back(stod(value));
            }
            catch (std::invalid_argument& e)
            {}
        }

    }

    return v;


}

Input Data:: resolveInput (std::string input)
{

    if( input == "GEN" ) return Generation;
    if( input == "REP" ) return Replicas;
    if( input == "SITES")return Sites;
    if( input == "MUT" ) return Mutations;

    return NoInput ;

}
void Data::setMutations(std::vector<double> list)
{
    for(auto mig  : list )
    {
        mutations.push_back(mig);
    }
}

std::vector<double> Data::getAllelesFq() const {
	return alleleFq;
}
