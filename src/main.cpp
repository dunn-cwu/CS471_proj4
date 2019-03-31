#include <iostream>
#include "proj1.h"
#include "inireader.h"

using namespace std;

int main(int argc, char** argv) 
{
   if (argc <= 1)
   {
      cout << "Error: Missing command line parameter." << endl;
      cout << "Proper usage: " << argv[0] << " [param file]" << endl;
      return EXIT_FAILURE;
   }

   proj1::mfuncExperiment ex;

   cout << "Input parameters file: " << argv[1] << endl;
   cout << "Initializing experiment ..." << endl;

   if (!ex.init(argv[1]))
      return EXIT_FAILURE;

   return ex.runAllFunc();
}
