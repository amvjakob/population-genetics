#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <vector>

/*!
  This is a random number class based on standard c++-11 generators
 */
class RandomDist {
    
public:
	/*!
	  Initializes the generator \ref rng with the Mersenne twister *mt19937* engine. 
	  Must provide mean *m*, standard deviation *sd* and sample size *ns*, 
	  set whether to use normal distribution (default uniform) and can provide a seed *s*.
	 */
    RandomDist(double m, double sd, int ns, bool n = false);
    
	/*!
	  Returns a vector of random doubles corresponding to the parameters set in the constructor.
	*/
    std::vector<double> generate_numbers();

	/** \brief Get a number following a binomial distribution
	 *
	 * Uses a mt19937 mersenne twister.
	 *
	 * */
    static int binomial(int n, double p);
    
    static int uniformIntSingle(int min, int max);
    static double uniformDoubleSingle(double min, double max);
    
    static void uniformIntVector(std::vector<int>& toFill, int min, int max);
    static void uniformDoubleVector(std::vector<double>& toFill, double min, double max);
    
    static void multinomial(std::vector<unsigned int>& pop, int n);
    static std::vector<unsigned int> multinomialByValue(const std::vector<unsigned int>& pop, int n);
     
private:
	static std::random_device rd;
	static std::mt19937 rng;

    void uniform(std::vector< double >&);
    void normal(std::vector< double >&);
    
    double mean, sd;
    int nsample;
    bool normdist;

};

#endif
