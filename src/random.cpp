#include <iostream>
#include "random.hpp"
#include <algorithm>

Randomdist::Randomdist(double m, double s, int ns, bool n, unsigned long int seed) 
    : mean(m), sd(s), nsample(ns), normdist(n) {
    if (s <= 0) {
        throw(1);
    }
    if (ns <= 0) {
        throw(2);
    } 
    if (seed == 0) {
        std::random_device rd;
        seed = rd();
    }
    rng = std::mt19937(seed);
}

std::vector< double > Randomdist::generate_numbers() {
    std::vector< double > result(nsample);
    if (normdist) {
        normal(result);
    }
    else {
        uniform(result);
    }
    return result;
}

void Randomdist::uniform(std::vector< double > &res) {
    
    double delta = sd*sqrt(3.0);
    double lower = mean-delta, upper = mean+delta;
    std::uniform_real_distribution<> unif(lower, upper);
    
    for ( auto I = res.begin(); I != res.end(); I++ ) {
        *I = unif(rng);
    }
}

void Randomdist::normal(std::vector< double > &res) {
   
    std::normal_distribution<> norm(mean, sd);
   
    for (auto I = res.begin(); I != res.end(); I++) {
        *I = norm(rng);
    }
}
