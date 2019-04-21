/**
 * @file main.cpp
 * @author Andrew Dunn (Andrew.Dunn@cwu.edu)
 * @brief Program entry point. 
 * Creates and runs CS471 project 1 experiment.
 * @version 0.1
 * @date 2019-04-01
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include <iostream>
#include <sstream>
#include "experiment.h"

using namespace std;

template<class T>
int runExp(const char* paramFile)
{
   // Create an instance of the project 1 experiment class
   mfunc::Experiment<T> ex;

   cout << "Datatype size: " << sizeof(T) << endl;
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

   int dataType = 1;

   if (argc > 2)
   {
      std::stringstream ss(argv[2]);
      ss >> dataType;
      if (!ss) dataType = 1;
   }

   if (dataType < 0 || dataType > 2)
   {
      cout << dataType << " is not a valid data type index. Value must be between 0 and 2." << endl;
      dataType = 1;
   }

   // Run experiment and return success code
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
