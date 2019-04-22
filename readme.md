Course: CS471
Assignment: Project 2
Student: Andrew Dunn
Instructor: Dr. Donald Davendra

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
in the [project]/source/params/ directory. To run Blind Search for
10, 20, and 30 dimensions with 30 iterations, open a terminal and execute 
the following commands:

```
cd [Path-to-/source-dir]
./unix-run-blindsearch.sh
```

Where [Path-to-/source-dir] is the file system path to the [project]/source/ directory
within this project. You should see a message that says all tests were ran after they complete.
Results files will be placed in the [project]/source/results directory.

--------

To run Local Search with 10, 20, and 30 dimensions with 30 iterations, open a terminal and 
execute the following commands:

```
cd [Path-to-/source-dir]
./unix-run-localsearch.sh
```
Where [Path-to-/source-dir] is the file system path to the [project]/source/ directory
within this project. You should see a message that says all tests were ran after they complete.
Results files will be placed in the [project]/source/results directory.

--------

To run both Blind Search and Local Search in 10, 20, and 30 dimensions with 30 iterations,
open a terminal and execute the following commands:

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
which will run the search algorithms with different dimensions. For example:

```
cd [Path-to-/source-dir]
./build/release/cs471_proj1.out ./param/param-blindsearch-10dim.ini
```

Available files:

./param/param-blindsearch-10dim.ini
./param/param-blindsearch-20dim.ini
./param/param-blindsearch-30dim.ini

./param/param-localsearch-10dim.ini
./param/param-localsearch-20dim.ini
./param/param-localsearch-30dim.ini

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
file which will run both Blind Search and Local Search for 30 iterations
in 10, 20, and 30 dimensions. Results files will be placed in the 
[project]/source/results directory.

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

The 'alpha' entry is used by Local Search when calculating the neighbor vector.

The 'algorithm' entry allows you to select which search algorithm to run.
0 = Blind Search and 1 = Local Search.

The 'results_file' entry is the file path (without spaces) to the file you 
wish to export the search algorithm fitness results to.

The 'exec_times_file' entry is the file path (without spaces) to the file you
wish to export the search algorithm execution times to.

--

The 'function_range' section lets you specify the random number
generator bounds for each function's data population.

---------------------------------
