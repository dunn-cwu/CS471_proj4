/**
 * @file main.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Program entry point. 
 * Creates and runs CS471 project 3 experiment.
 * @version 0.2
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include <iostream>
#include <sstream>
#include "experiment.h"

using namespace std;

/**
 * @brief Runs the experiment using the given data type
 * and parameter file. Currently supports three different
 * data types: float, double, and long double.
 * 
 * @tparam T 
 * @param paramFile 
 * @return int 
 */
template<class T>
int runExp(const char* paramFile)
{
   // Create an instance of the project 1 experiment class
   mfunc::Experiment<T> ex;

   // Print size of selected data type in bits
   cout << "Float size: " << (sizeof(T) * 8) << "-bits" << endl;
   cout << "Input parameters file: " << paramFile << endl;
   cout << "Initializing experiment ..." << endl;

   // If experiment initialization fails, return failure
   if (!ex.init(paramFile))
      return EXIT_FAILURE;
   else
      return ex.testAllFunc();
}

int main(int argc, char** argv) 
{
   // Make sure we have enough command line args
   if (argc <= 1)
   {
      cout << "Error: Missing command line parameter." << endl;
      cout << "Proper usage: " << argv[0] << " [param file]" << endl;
      return EXIT_FAILURE;
   }

   // Default data type is double
   int dataType = 1;

   // User specified a data type, retrieve the value
   if (argc > 2)
   {
      std::stringstream ss(argv[2]);
      ss >> dataType;
      if (!ss) dataType = 1;
   }

   // Verify specified data type switch
   if (dataType < 0 || dataType > 2)
   {
      cout << dataType << " is not a valid data type index. Value must be between 0 and 2." << endl;
      dataType = 1;
   }

   // Run experiment with correct data type and return success code
   switch (dataType)
   {
      case 0:
         return runExp<float>(argv[1]);
      case 1:
         return runExp<double>(argv[1]);
      case 2:
         return runExp<long double>(argv[1]);
      default:
         return EXIT_FAILURE;
   }
}

// =========================
// End of main.cpp
// =========================
