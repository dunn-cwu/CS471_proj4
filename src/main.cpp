#include <iostream>
#include "proj1.h"

int main(int argc, char *argv[]) {
   std::vector<double> resultsV;
   double time;

   proj1::mfuncExperiment ex(10, 30);
   ex.runFunc(1, resultsV, time);
   std::cout << time << std::endl;

   return 0;
}
