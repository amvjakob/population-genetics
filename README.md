# Team 4 Genetics Project

### Aim
Simulate population genetics.

The first version of the project is to simulate a basic version of the
Wright Fisher model. This consists of having random alleles that have the same
probability of integrating the next generation. The simulation is a model of the
Genetic drift within a population.

An additional feature of the simulation is to have a more selective allele repartition
and a consideration of additional factors (migration, mutations, selection) that
will create a more realistic simulation of population genetics.

As such, the program can be launched by using various modes. Each mode corresponds to a
specific model that will take into consideration one of the following factors :
migration, selection, mutation and time dependent population size (bottleneck
effect). Note that the mutation mode requires additional precision from the user
to decide on a specific mutation model (Cantor, Felsenstein, Kimura).

The program can be launched using two different methods:
  1- using a fasta file
  2- using the input.txt file were the different parameters (population size,
    number of generations, mutation marker sites, etc.) can be modified for the
    desired simulation

Finally, the results.txt file generated for each simulation can be used to generate
various graphs using jupyterNotebook and plotly library.
