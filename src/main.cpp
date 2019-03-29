#include <iostream>
#include "mfunc.h"

int main(int argc, char *argv[]) {
   double testV[5] = { 5.4, 67.2, 26.2, -123.6, 400.3 };

   double result = mfunc::schwefel(testV, 5);

   std::cout << result << std::endl;

   return 0;
}
