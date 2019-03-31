#include <iostream>
#include "proj1.h"

int main(int argc, char *argv[]) {
   proj1::mfuncExperiment ex(30, 30);
   ex.runAllFunc("data.csv");
   return 0;
}
