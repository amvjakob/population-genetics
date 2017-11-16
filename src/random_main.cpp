#include "random.hpp"
#include <iostream>
#include <tclap/CmdLine.h>

Randomdist* parse_args(int argc, char **argv) {
    TCLAP::CmdLine cmd("Random number generator");
    TCLAP::ValueArg<int> nsample("N", "sample_size", "Number of random numbers to generate",  true, 1, "int");
    cmd.add(nsample);
    TCLAP::ValueArg<double> meanval("m", "mean", "Mean value",  true, 0, "double");
    cmd.add( meanval );
    TCLAP::ValueArg< double > sdev("s", "std_dev", "Standard deviation",  true, 1, "double");
    cmd.add( sdev );
    TCLAP::SwitchArg  uniform("u", "uniform", "Uniform distribution", true);
    cmd.add( uniform );
    TCLAP::SwitchArg  normal("n", "normal", "Normal distribution", true);
    cmd.add(normal);
    cmd.parse( argc, argv );
    double m = meanval.getValue(), s = sdev.getValue();
    int ns = nsample.getValue();
    bool norm =  normal.getValue();
    if ( !(norm ^ uniform.getValue()) ) {
        throw(3);
    }
    Randomdist* rng = new Randomdist(m, s, ns, norm);
    return rng;
}

int main(int argc, char **argv) {
    Randomdist *rng = 0;
    try {
        rng = parse_args(argc, argv);
    } catch(const int n) {
        switch (n) {
        case 1:
            std::cerr << "Standard deviation must be positive\n";
            break;
        case 2:
            std::cerr << "Sample size must be positive\n";
            break;
        case 3:
            std::cerr << "Please choose either normal or uniform distribution\n";
            break;
        }
        return n;
    }
    for (auto I : rng->generate_numbers()) {
        std::cout << I << " "; std::cout << "\n";
    }
    if (rng) delete rng;
    return 0;
}

