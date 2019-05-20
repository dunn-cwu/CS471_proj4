Course: CS471
Assignment: Project 4
Student: Andrew Dunn
Instructor: Dr. Donald Davendra

---------------------------------
Project Description
---------------------------------

This project tests 18 different dimension scalable 
mathematical functions with a random population data
set. The goal is to find a population that produces
a fitness close to the optimal value using various
search algorithms. The search algorithms used are
particle swarm, firefly algorithm, and harmony search.

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
in the [project]/source/params/ directory. To run particle swarm for
example, open a terminal and execute the following commands:

```
cd [Path-to-/source-dir]
./unix-run-pso.sh
```

Where [Path-to-/source-dir] is the file system path to the [project]/source/ directory
within this project. You should see a message that says all tests were ran after they complete.
Results files will be placed in the [project]/source/results directory.

--------

To run the firefly or harmony search algorithms, open a terminal and 
execute the following commands:

```
cd [Path-to-/source-dir]
./unix-run-fa.sh
./unix-run-hs.sh
```
Where [Path-to-/source-dir] is the file system path to the [project]/source/ directory
within this project. You should see a message that says all tests were ran after they complete.
Results files will be placed in the [project]/source/results directory.

--------

To run both all three search algorithms,
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
./build/release/cs471-proj4.out [Input-parameter-file]
```

Where [Input-parameter-file] is the path to the input parameter ini file.
There are three already provided within the [project]/source/params/ directory 
which will run the search algorithms with different parameters. For example:

```
cd [Path-to-/source-dir]
./build/release/cs471-proj4.out ./params/fa-1.ini
```

Available files:

Used for debugging:
./params/debug_params.ini

Firefly algorithm, four parameter strategies:
./params/fa-1.ini
./params/fa-2.ini
./params/fa-3.ini
./params/fa-4.ini

Harmony search, four parameter strategies:
./params/hs-1.ini
./params/hs-2.ini
./params/hs-3.ini
./params/hs-4.ini

Particle swarm, four parameter strategies
./params/pso-1.ini
./params/pso-2.ini
./params/pso-3.ini
./params/pso-4.ini

---------------------------------
Optional Run Command Line Argument
---------------------------------

When running the project manually, you can specify an additional command line
argument that allows you to select a specific data type to be used in the experiment.
To do this, open a terminal and execute the following commands:

```
./build/release/cs471-proj4.out [Input-parameter-file] [Datatype-id]
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
It contains two different main sections, 'test' and 'function_range'.

The 'test' section contains various settings to control how the
experiment is ran and which files are produced. The 'function_range'
section lets you specify the random number generator bounds for each 
function's data population.

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
0 = Particle swarm, 1 = Firefly algorithm, and 2 = Harmony search.

The 'results_file' entry is the file path (without spaces) to the file you 
wish to export the best fitness results to.

The 'worst_fit_file' entry is the file path (without spaces) to the file you 
wish to export the worst fitness results to.

The 'exec_times_file' entry is the file path (without spaces) to the file you
wish to export the search algorithm execution times to.

The 'func_calls_file' entry is the file path (without spaces) to the file you
wish to export objective function call counts to.

--

The 'particle_swarm' section lets you specify parameters specific to the particle swarm algorithm:

The 'c1' entry specifies c1 random chance for particle swarm (0-1.0)

The 'c2' entry specifies c2 random chance for particle swarm (0-1.0)

The 'k' entry specifies the velocity dampening factor for particle swarm

--

The 'firefly' section lets you specify parameters specific to the firefly algorithm:

The 'alpha' entry specifies the alpha parameter for the firefly algorithm

The 'betamin' entry specifies the betamin parameter for the firefly algorithm

The 'gamma' entry specifies the gamma parameter for the firefly algorithm

--

The 'harmony_search' section lets you specify parameters specific to the harmony search algorithm:

The 'hmcr' entry specifies the Harmony Memory Considering Rate parameter for the harmony search algorithm

The 'par' entry specifies the Pitch Adjusting Rate parameter for the harmony search algorithm

The 'bw' entry specifies the Bandwidth parameter for the harmony search algorithm

---------------------------------
