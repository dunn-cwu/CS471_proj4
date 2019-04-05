Course: CS471
Assignment: Project 1
Student: Andrew Dunn
Instructor: Dr. Donald Davendra

---------------------------------
Project Description
---------------------------------

This project tests 18 different dimension scalable 
mathematical functions with a random population data
set. The goal is to find a population that produces
a fitness close to the optimal value.

---------------------------------
Build Instructions
---------------------------------

Project Requirements:

A computer running Linux or MacOS
Cmake version 3.1.0+
g++ on Linux, or clang on MacOS

Your C++ compiler must support the C++11 or greater standard.

----

Within the source/ directory you will find a few different
shell scripts. To build this project you can run either
'build-debug.sh' or 'build-release.sh' from within the
source/ folder. To build the release binary, open a terminal
and execute the following commands:

cd [Path-to-/source-dir]
./build-release.sh

Where [Path-to-/source-dir] is the file system path to the /source directory
within this project. This script will run cmake and make, building the 
project with your default C++ compiler. The produced binary will be moved to
the build/release directory.

---------------------------------
Run Instructions
---------------------------------

To run this project, a simple shell script is provided which will
run the release binary three times with the three different input
parameter files included in the source directory. To run this script,
open a terminal and execute the following commands:

cd [Path-to-/source-dir]
./run-all.sh

Where [Path-to-/source-dir] is the file system path to the /source directory
within this project. You should see a message that says all tests were ran.

----

If you wish to run this project manually, within the terminal run the
following command:

/build/release/cs471_proj1.out [Input-parameter-file]

Where [Input-parameter-file] is the path to the input parameter ini file.
There are three already provided within the /source directory which
will run the experiment with different dimensions. These files are:

iparam-10dim.ini
iparam-20dim.ini
iparam-30dim.ini

---------------------------------
Input parameter file format
---------------------------------

The input parameter file is a configuration file in the *.ini format.
It contains two different sections, 'test' and 'function_range'.

The 'test' section contains various settings to control how the
experiment is ran and which files are produced.

--

Within the 'test' section:

The 'population' entry sets the size of the population, i.e. how
many times each funcion is executed.

The 'dimensions' entry sets the number of dimensions for each of the
population vectors that are passed to the functions.

If you wish all population data to be output to files, set 
'output_population' to true.

If you wish all resulting fitness data to be output to files, set
'output_fitness' to true.

'results_file' is the file path (without spaces) to the file you 
wish to contain the experiment results and analysis table.

--

The 'function_range' section lets you specify the random number
generator bounds for each function's data population.

---------------------------------

