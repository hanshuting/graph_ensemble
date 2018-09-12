#include <iostream>

#include "common.hpp"
#include "info.hpp"
#include "configuration.hpp"
#include "veclabelvecfeaturedc.hpp"

using namespace std;

int main(int argc, char** argv)
{
   CInfo info;
   Configuration &config = Configuration::GetInstance();

   try {
      config.ReadFromFile(argv[1]);
      
      CVecLabelVecFeatureDC ds(&info);
      
      Scalar norm1 = 0, norm2 = 0, norminf = 0;
      
      ds.X->Norm1(norm1);
      ds.X->Norm2(norm2);
      ds.X->NormInf(norminf);
      cout << "L_1 norm of feature matrix:       " << norm1 << endl;
      cout << "Frobenius norm of feature matrix: " << norm2 << endl;
      cout << "L_\\infty norm of feature matrix:  " << norminf << endl;
   }
   catch(CBMRMException e) {
      cout << e.Report() << endl;
   }

   return EXIT_SUCCESS;
}
