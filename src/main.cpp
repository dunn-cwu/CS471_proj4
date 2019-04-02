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
#include "proj1.h"

using namespace std;

int main(int argc, char** argv) 
{
   // Make sure we have enough command line args
   if (argc <= 1)
   {
      cout << "Error: Missing command line parameter." << endl;
      cout << "Proper usage: " << argv[0] << " [param file]" << endl;
      return EXIT_FAILURE;
   }

   // Create an instance of the project 1 experiment class
   proj1::mfuncExperiment ex;

   cout << "Input parameters file: " << argv[1] << endl;
   cout << "Initializing experiment ..." << endl;

   // If experiment initialization fails, return failure
   if (!ex.init(argv[1]))
      return EXIT_FAILURE;

   // Run experiment and return success code
   return ex.runAllFunc();
}

// =========================
// End of main.cpp
// =========================
