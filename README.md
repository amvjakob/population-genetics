# Genetics Project- Team 4 (EPFL)

## Aim
Simulate population genetics.

The program developed aims to simulate the Wright Fisher model. This consists in having random alleles that have the same probability of integrating the next generation. The simulation is a model of the genetic drift within a population.

An additional feature of the simulation is to have a more selective allele distribution to create a more realistic simulation of population genetics, by taking into consideration specific factors such as:
1. Mutations
2. Migrations
3. Selection
4. Time-dependent population

As such, the program can be launched by using various modes. Each mode corresponds to a specific model that will take into consideration one of the aforementioned parameters.
Note that the mutation mode requires additional precision from the user to decide on a specific mutation model (Jukes-Cantor, Kimura, Felsenstein) and the
migration mode requires the pattern followed by the sub-groups (complete graph, ring, star) as well as the allele repartition (random or input by the user).

## Setup the project
To run the program, follow these steps:
1. Download or clone the repository
2. Navigate to the cloned / downloaded folder (Team_4)
3. `mkdir build`, `cd build` to enter the build folder
4. `cmake ..` to run CMake and generate the makefiles
5. `make` to make both the simulation as well as the tests.
6. `make test` to make and run the tests
7. `make doc` to generate the Doxygen documentation


The program can be launched using two different methods:
1. `./Genetics`, runs the simulation using only the default input file (data/input.txt) with the desired parameters for the simulation
2. `./Genetics path/to/input.txt`, runs the simulation with the specified input file
3. `./Genetics path/to/input.txt path/to/fasta.fa`, runs the simulation with the specified input file and the given fasta file. The input file is modified by the user to decide on the various simulation parameters (population size, number of generations, mutation marker sites, etc.).
 
Finally, the `results.txt` file generated for each simulation and can be used to generate various graphs using jupyterNotebook and plotly library.
 
#### Note
For a simulation using mutation models, a fasta file is mandatory.


## Special Feature:
The program is coded using multiple threads. Each thread executes a single simulation, allowing for replicas to run simultaneously. This allows for faster simulation.

## Run a simulation:
The user must use the input.txt file to choose:
- the number of generations
- the population size
- the execution mode (mutation, migration, selection, time-dependent population)
- the marker sites (if using the mutation model)

If working with the mutation model, the user must precise:
- the mutation model with the appropriate probabilities

If working with the migration mode, the user must precise:
- the mutation pattern
- the migration rates

If working with the selection model, the user must precise:
- the selection rates for each allele

If working with the time-dependent population model, the user must precise:
- the time (simulation step) at which the population change occurs
- the time (simulation step) at which it ends
- the factor by which the population will change (a factor of 2 means the population will be halved during the specified time interval)

## Authors
- Loïc Bruchez
- Sofia Dandjee
- Aurélie Ducrot
- Anthony Jakob
- Kamil Seghrouchni
- Luca Zanetti

## Presentation
Presentation of the project (PPTX): <https://prezi.com/p/jfcpou3kqbcd/>