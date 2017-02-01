#ifndef _INFO_HPP_
#define _INFO_HPP_

#include "common.hpp"
#include <string>

using namespace std;

/**   Class to keep program-wide information such as parameters, precision tolerance, etc.
 */
namespace INFO {
   enum PROGRAM_MODE {LEARNING, PREDICTION};
   enum COMPUTATION_MODE {SERIAL, PARALLEL_CENTRALISED_DATA, PARALLEL_DISTRIBUTED_DATA};
}

class CInfo {
   public:
      int    verbosity;             // verbosity of the program
      int    convergenceLog;        // keep convergence statistics in files?
      int    programMode;           // learning or predicting?
      int    lossFunctionType;      // the type (#defined number) of the risk minimization problem
      int    regularizerType;       // the type of solver/regularizer e.g. L1, L2, ...
      int    bias;                  // with additional bias feature (append as new constant dimension of data)
      int    dimW;                  // dimensionality weight vector (and data)
      int    numOfHyperplane;       // number of hyperplane: in binary classification and regression, it's 1; in n-class multiclass classification, it's n.
      double lambda;                // regularization constant
      double C;                     // margin-error trade off constant
      double epsilonTol;            // termination tolerance for epsilon termination criterion (SVMStruct)
      double gammaTol;              // termination tolerance for gamma termination criterion (BMRM)
      int    maxNumOfIter;          // maximum number of iteration allowed 
      
      string hostname;              // name of the machine running this program
      string datasetFn;             // dataset filename (the feature vector part)
      string labelsetFn;            // labelset filename (the label part)
      string datasetList;           // file containing a list of sub-dataset filenames
      string modelFn;               // file for storing model (weight vector and other information such as dimX, numOfHyperplane, etc.)
      
      // for prediction phase
      int    outputType;            // type of prediction output: FVal only or label only or both
      string outputLabelFn;         // predicted label output filename
      string outputFValFn;          // decision function value filename

      // parallel and distributed computation related
      int    procID;                // ID of machine running this process
      int    numOfNode;             // number of machines running bmrm:
                                    //   1 := serial bmrm; >1 := parallel bmrm
      int    numOfUsedNode;         // in the cases where numOfNode > nSample, only a portion of the 
                                    //   available machines will be used.   
      int    computationMode;       // computation mode: {serial, parallel/distributed}
                                    //   1 := serial (with centralised dataset)
                                    //   2 := parallel with centralised dataset
                                    //   3 := parallel with distributed dataset  
      
      // constructors
      CInfo();
      CInfo(int argc, char** argv) {CInfo();}  // TODO: implement this!

      // destructor
      virtual ~CInfo() {}

      // method
      void Print();
};

#endif
