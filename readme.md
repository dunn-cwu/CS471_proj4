Course: CS471
Assignment: Project 3
Student: Andrew Dunn
Instructor: Dr. Donald Davendra

---------------------------------
Important note from student
---------------------------------

I had a difficult time with this project. I ran into
multiple roadblocks when implementing the two algorithms,
especially the genetic algorithm which set me back a few days.
I was getting wierd results that took me a long time to figure out
what I was doing wrong, and I think my DE implementation is still
not working 100% as it should.

In order to leave enough time for the sizable report, I decided to 
leave out the additional GA selection algorithms and multiple 
crossover function to focus on a working code. I apologize for my failings, 
but it is what it is. At the end of the day I had to make a call in order 
to finish the project on time, and so this was my choice.

---------------------------------
Project Description
---------------------------------

This project tests 18 different dimension scalable 
mathematical functions with a random population data
set. The goal is to find a population that produces
a fitness close to the optimal value using various
search algorithms.

---------------------------------
Project Requirements
---------------------------------

1. A computer running Linux, MacOS, or Windows.
2. Some version of g++ that supports C++11 (g++, MinGW, or clang)
3. Cmake 3.1+

----

This program was succesfully built, ran, and tested on the following
operating systems and tool chains:

Ubuntu 19.10
g++ version 8.2 and cmake 3.10

MacOS 10.13
clang 10.0.0 and cmake 3.14

Windows 10
g++ via MinGW-W64 8.1 and cmake 3.14

Your C++ compiler must support the C++11 standard or greater.

---------------------------------
Build Instructions - Unix Based Machines
---------------------------------

For Ubuntu and MacOS:

Within the [project]/source/ directory you will find a few different
shell scripts. To build this project you can run either
'unix-build-debug.sh' or 'unix-build-release.sh' from within the
[project]/source/ directory. To build the release binary, open a terminal
and execute the following commands:

```
cd [Path-to-/source-dir]
./build-release.sh
```

Where [Path-to-/source-dir] is the path to the [project]/source/ directory
within this project. This script will run cmake and make, building the 
project with your default C++ compiler which should be auto-detected. 
The produced binary will be moved to the build/release directory.

---------------------------------
Build Instructions - Windows Based Machines
---------------------------------

For Windows 10 and earlier:

Within the [project]\source\ directory you will find a few different
batch scripts. To build this project in release mode you can run win-build-release.bat 
by simply double clicking it. This script will run MinGW cmake and make, building the 
project with your default MinGW C++ compiler which should be auto-detected if your
PATH system variable is set up correctly. The produced binary will be moved to the 
build/release directory.

---------------------------------
Run Instructions - Unix Based Machines
---------------------------------

For Ubuntu and MacOS:

To run this project, a few shell scripts are provided that will
run the release binary using different test parameter files included 
in the [project]/source/params/ directory. To run genetic algorithm for
50 iterations with a population size of 200 in 30 dimenstions, open a
terminal and execute the following commands:

```
cd [Path-to-/source-dir]
./unix-run-GA.sh
```

Where [Path-to-/source-dir] is the file system path to the [project]/source/ directory
within this project. You should see a message that says all tests were ran after they complete.
Results files will be placed in the [project]/source/results directory.

--------

To run differential evolution for all 10 strategies with population sizes of 200 in 30 dimensions, open a terminal and 
execute the following commands:

```
cd [Path-to-/source-dir]
./unix-run-DE.sh
```
Where [Path-to-/source-dir] is the file system path to the [project]/source/ directory
within this project. You should see a message that says all tests were ran after they complete.
Results files will be placed in the [project]/source/results directory.

--------

To run both genetic algorithm and differential evolution with the parameters specified above,
open a terminal and execute the following command:

```
cd [Path-to-/source-dir]
./unix-run-all.sh
```
Where [Path-to-/source-dir] is the file system path to the [project]/source/ directory
within this project. You should see a message that says all tests were ran after they complete.
Results files will be placed in the [project]/source/results directory.

--------

If you wish to run this project manually, within the terminal run the
following command:

```
./build/release/cs471_proj1.out [Input-parameter-file]
```

Where [Input-parameter-file] is the path to the input parameter ini file.
There are three already provided within the [project]/source/params/ directory 
which will run the search algorithms with different parameters. For example:

```
cd [Path-to-/source-dir]
./build/release/cs471_proj1.out ./params/GA_roulette.ini
```

Available files:

./params/GA_roulette.ini

./params/DE_Strat1.ini
./params/DE_Strat2.ini
./params/DE_Strat3.ini
./params/DE_Strat4.ini
./params/DE_Strat5.ini
./params/DE_Strat6.ini
./params/DE_Strat7.ini
./params/DE_Strat8.ini
./params/DE_Strat9.ini
./params/DE_Strat10.ini

---------------------------------
Optional Run Command Line Argument
---------------------------------

When running the project manually, you can specify an additional command line
argument that allows you to select a specific data type to be used in the experiment.
To do this, open a terminal and execute the following commands:

```
./build/release/cs471_proj1.out [Input-parameter-file] [Datatype-id]
```

Where [Datatype-id] is an integer ranging from 0 to 2:

0 = 32-bit float
1 = 64-bit float
2 = 128-bit float (when supported by your system)

---------------------------------
Run Instructions - Windows based machines
---------------------------------

A few batch scripts are provided in the [project]\source\ directory.

To run this project in Windows, simply double click the the win-run-all.bat
file which will run both the genetic algorithm and differential evolutionary algorithms 
tests. Results files will be placed in the [project]/source/results directory.

---------------------------------
Input parameter file format
---------------------------------

The input parameter file is a configuration file in the *.ini format.
It contains two different sections, 'test' and 'function_range'.

The 'test' section contains various settings to control how the
experiment is ran and which files are produced.

--

Within the 'test' section:

The 'population' entry sets the number of rows in the population vector matix.

The 'dimensions' entry sets the number of dimensions for each of the
population vectors.

The 'iterations' entry sets the number of test iterations for the selected
search algorithm.

The 'num_threads' entry sets the number of worker threads you want to use
to run the experiment. Note that you want to set this value to be equal or
close to the number of CPU's/CPU cores available in your system.

The 'algorithm' entry allows you to select which search algorithm to run.
0 = Genetic Algorithm and 1 = Differential Evolution.

The 'results_file' entry is the file path (without spaces) to the file you 
wish to export the search algorithm fitness results to.

The 'exec_times_file' entry is the file path (without spaces) to the file you
wish to export the search algorithm execution times to.

--

The 'genetic_alg' section lets you specify parameters specific to the genetic algorithm:

The 'generations' entry specifies the number of generations to run the algorithm for.

The 'crossover_prob' entry specifies the crossover probability (0-1.0)

The 'mutation_prob' entry specifies the mutation probability (0-1.0)

The 'mutation_range' entry specifies the mutation range

The 'mutation_precision' entry specifies the mutation precision

The 'elitism_rate' entry specifies the elitisn rate (0-1.0)

--

The 'differential_evo' section lets you specify parameters specific to the differential evolution algorithm:

The 'generations' entry specifies the number of generations to run the algorithm for.

The 'crossover_prob' entry specifies the crossover probability (0-1.0)

The 'scalefactor_1' entry specifies the primary mutation scale factor

The 'scalefactor_2' entry specifies the secondary (lambda) mutation scale factor

The 'strategy' entry specifies the mutation and crossover strategy to be used. There are 10 different strategies that
can be selected using an index 0-9. These strategies are:

0 = Best1Exp
1 = Rand1Exp
2 = RandToBest1Exp
3 = Best2Exp
4 = Rand2Exp
5 = Best1Bin
6 = Rand1Bin
7 = RandToBest1Bin
8 = Best2Bin
9 = Rand2Bin

--

The 'function_range' section lets you specify the random number
generator bounds for each function's data population.

---------------------------------
