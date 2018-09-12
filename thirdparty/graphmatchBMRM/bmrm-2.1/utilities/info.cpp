#ifndef _INFO_CPP_
#define _INFO_CPP_

#include "common.hpp"
#include "info.hpp"

CInfo::CInfo()
   : verbosity(0),
     convergenceLog(0),
     programMode(INFO::LEARNING),
     lossFunctionType(0),
     regularizerType(0),
     bias(0),
     dimW(0),
     numOfHyperplane(1),
     lambda(1),
     C(1),
     epsilonTol(1e-3),
     gammaTol(-1),
     maxNumOfIter(2000),
     hostname(""),
     datasetFn(""),
     labelsetFn(""),
     datasetList(""),
     modelFn(""),
     outputType(LABEL_AND_FVAL),
     outputLabelFn(""),
     outputFValFn(""),
     procID(0),
     numOfNode(1),
     numOfUsedNode(1),
     computationMode(INFO::SERIAL)
{}


/** Print (selected) information
 */
void CInfo::Print()
{
  printf("=====  Information and parameters   =====\n");
  printf(" verbosity             : %d\n",    verbosity);
  printf(" convergenceLog        : %d\n",    convergenceLog);
  printf(" loss function type    : %d\n",    lossFunctionType);
  printf(" regularizer           : %d\n",    regularizerType);
  printf(" added bias feature?   : %d\n",    bias);
  printf(" num of hyperplane     : %d\n",    numOfHyperplane);
  printf(" dimW                  : %d\n",    dimW);
  printf(" lambda                : %.6e\n",  lambda);
  printf(" C                     : %.6e\n",  C);
  printf(" Machine name          : %s\n",    hostname.c_str());
  printf(" Data file             : %s\n",    datasetFn.c_str());
  printf(" Label file            : %s\n",    labelsetFn.c_str());
  printf(" Model file            : %s\n",    modelFn.c_str());
  printf(" Max. Col. Gen. Iter   : %d\n",    maxNumOfIter);
  
  double tol = epsilonTol;
  string g_or_e = "(epsilon)";
  if(gammaTol > 0) {
     tol = gammaTol;
     g_or_e     = "(gamma)  ";
  }
  printf(" Tolerance %s : %.8f\n",           g_or_e.c_str(), tol);
  printf(" # of avail. machines  : %d\n",    numOfNode);
  printf(" # of machines used    : %d\n",    numOfUsedNode);
  printf(" Computation mode      : %d\n",    computationMode);
  printf("=========================================\n");
}
   
#endif
