#include "Random.hpp"
#include <algorithm>

std::random_device RandomDist::rd;
std::mt19937 RandomDist::rng = std::mt19937(RandomDist::rd());


RandomDist::RandomDist(double m, double s, int ns, bool n) 
    : mean(m), sd(s), nsample(ns), normdist(n) {
    if (s <= 0) {
        throw(1);
    }
    if (ns <= 0) {
        throw(2);
    }
}

std::vector<double> RandomDist::generate_numbers() {
    std::vector<double> result(nsample);
    if (normdist) {
        normal(result);
    }
    else {
        uniform(result);
    }
    return result;
}

int RandomDist::binomial(int n, double p) {
	std::binomial_distribution<int> dbinom(n, p);
	return dbinom(rng);
}

int RandomDist::uniformIntSingle(int min, int max) {
	// init random distribution
	std::uniform_int_distribution<int> distr(min, max);
	
	// get one value
	return distr(RandomDist::rng);
}

double RandomDist::uniformDoubleSingle(double min, double max) {
	// init random distribution
	std::uniform_real_distribution<> distr(min, max);
	
	// get one value
	return distr(RandomDist::rng);
}

void RandomDist::uniformIntVector(std::vector<int>& toFill, int min, int max) {
	// init random distribution
	std::uniform_int_distribution<int> distr(min, max);

	// fill table with generated values
	std::generate(
		toFill.begin(),
		toFill.end(), 
		[&]() { 				
			return distr(RandomDist::rng); 
		}
	); 
}

void RandomDist::uniformDoubleVector(std::vector<double>& toFill, double min, double max) {
	// init random distribution
	std::uniform_real_distribution<> distr(min, max);

	// fill table with generated values
	std::generate(
		toFill.begin(),
		toFill.end(), 
		[&]() { 				
			return distr(RandomDist::rng); 
		}
	); 
}

void RandomDist::uniform(std::vector< double > &res) {
    
    double delta = sd * sqrt(3.0);
    double lower = mean - delta, upper = mean + delta;
    std::uniform_real_distribution<> unif(lower, upper);
    
    for (auto I = res.begin(); I != res.end(); I++) {
        *I = unif(rng);
    }
}

void RandomDist::normal(std::vector< double > &res) {
   
    std::normal_distribution<> norm(mean, sd);
   
    for (auto I = res.begin(); I != res.end(); I++) {
        *I = norm(rng);
    }
}
