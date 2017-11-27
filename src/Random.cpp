#include "Random.hpp"
#include <algorithm>
#include <cassert>

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



void RandomDist::multinomial(std::vector<unsigned int>& pop, int n) {
	int total = 0;
	for (auto& count : pop)
		total += count;
	
	for (auto& count : pop) {
		// remaining parent population size should be 0 or more	
		assert(total >= 0);
		
		double p;
		
		if (total == 0) {
			// the only way for the parent population to be 0 is if the 
			// current allele is not present anymore (wiped out)
			assert(total == 0);
			p = 0;
		} else {
			// generate new allele copy number
			p = count * 1.0 / total;
		}
		
		// reduce residual "gene pool"
		total -= count;
		
		// generate new number of allele copies in population
        int newCount = RandomDist::binomial(n, p);
		count = newCount;
		
		// reduce residual population size
		n -= newCount;
	}
	
	assert(n == 0);
	assert(total == 0);
}

/*
int RandomDist::hypergeometric(int m, int n, int k) {
	int N = m + n;
	int res = 0;
	
	assert(k <= N);
	
	std::uniform_int_distribution<int> distr;
	
	for (int i = 0; i < k; ++i) {
		distr = std::uniform_int_distribution<int>(0, N);
		
		if (distr(RandomDist::rng) < m) {
			++res;
			--m;
		}
		
		--N;
	}
	
	return res;
}

std::vector<int> RandomDist::multivariateGeometric(const std::vector<int>& population, int k) {
	// count population size
	int N = 0;
	for (auto& popCount : population)
		N += popCount;
		
	assert(N >= k);

	// count number of different types in population
	int m = (int) population.size();
	assert(m > 1);


	int nOther = N - population[0];
	
	std::vector<int> selected(m, 0);
	
	selected[0] = RandomDist::hypergeometric(population.at(0), nOther, k);
	
	for (int i = 1; i < m; ++i) {
		nOther -= population.at(i);
		
		k -= selected[i - 1];
		
		selected[i] = RandomDist::hypergeometric(population.at(i), nOther, k);
	}
	
	selected[m - 1] = k - selected[m - 2];
	
	return selected;
}
* */