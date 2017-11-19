#include "data.hpp"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cassert>

using namespace std;

Data::Data() :  dataName ("../doc/input.txt"),fastaName("../doc/test.fa"),populationSize(0), numberGenerations (0), replicates (0),numberAlleles(0) {

	assert(alleleFq.empty());
	assert(markerSites.empty());
	assert (sequences.empty());

	/*
	cout << "Please enter the name of your input file: " << endl;
	cin >> dataName;
	dataName += ".txt";

	cout << "Please enter the name of your fasta file: " << endl;
	cin >> fastaName;
	dataName += ".fa";
	*/
}



void Data::collectAll() {

	ifstream dataFile (dataName, ios::in);
	ifstream fastaFile (fastaName, ios::in);

	if (dataFile.is_open()) {

		collectDataFile(dataFile);
	} else {
		cout << "Ouverture  de l'input impossible" << endl;
	}

	if (fastaFile.is_open()) {

		collectFastaFile(fastaFile);

	} else {
		cout << "Ouverture  du fasta impossible" << endl;
	}

}

void Data::collectDataFile(ifstream& file) {


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


                case MigrationA :
                {
                    setMigrations(compar_son(2,"AD",line));
                }

                case MigrationT :
                {
                    setMigrations(compar_son(2,"TH",line));

                }

                case MigrationC :
                {
                    setMigrations(compar_son(2,"CY",line));
                }

                case MigrationG :
                {
                    setMigrations(compar_son(2,"GU",line));

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

std::list<double> Data::getMarkerSites() const {
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

list <double> Data::getMigrations() const {
	return migrations;
}

double Data:: comparison(double size,string s2,string l)
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




std::list<double> Data::compar_son(double size,std::string s2,std::string l)
{
    string key;
    string value ;
    std::list<double> v;

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
    if( input == "AD" ) return MigrationA;
    if( input == "TH" ) return MigrationT;
    if( input == "CY" ) return MigrationC;
    if( input == "GU" ) return MigrationG;

    return NoInput ;



}
void Data::setMigrations(std::list<double> list)
{
    for(auto mig  : list )
    {
        migrations.push_back(mig);
    }
}