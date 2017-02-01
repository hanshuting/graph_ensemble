#ifndef _COMMON_HPP_
#define _COMMON_HPP_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// #define's
// program mode
#define LEARNING_MODE      1    // do learning
#define PRED_AND_EVAL_MODE 2    // do prediction and evaluation (when labelset is supplied)
#define PREDICTION_MODE    3    // do prediction only (labelset is not supplied)

// regularizer type
#define L1NORM   1
#define L2NORMSQ 2

// BMRM learning type (loss functions)
#define SOFTMARGIN_LOSS            1
#define SQUARED_SOFTMARGIN_LOSS    2
#define HUBER_HINGE_LOSS           3
#define LOGISTIC_LOSS              4
#define EXPONENTIAL_LOSS           5
#define ROC_SCORE_LOSS             6
#define FBETA_SCORE_LOSS           7

#define MULTICLASSWTA_LOSS         8

#define ESVR_LOSS                  9

#define HUBER_LOSS                 10
#define LEAST_SQUARES_LOSS         11
#define LEAST_DEVIATION_LOSS       12

#define QUANTILE_LOSS              13
#define NOVELTY_LOSS               14
#define POISSON_LOSS               15

#define NDCG_LOSS                  16
#define MRR_LOSS                   17
#define WTA_LOSS                   18
#define CROSS_ENTROPY_LOSS         19
#define ORDINAL_REGRESSION_LOSS    20
#define PREFERENCE_RANKING_LOSS    21

#define NUMBER_OF_LOSS             21  // how many loss function we have


// BMRM prediction output type
#define NO_OUTPUT          0
#define LABEL_ONLY         1
#define FVAL_ONLY          2
#define LABEL_AND_FVAL     3

// computation mode
#define SERIAL_CENTRALISED_DS    1
#define PARALLEL_CENTRALISED_DS  2
#define PARALLEL_DISTRIBUTED_DS  3


// utility functions
bool IsBlankLine(std::string &line);

void trim(std::string& str);

void tokenize(const std::string& str, 
              std::vector<std::string>& tokens, 
              const std::string& delimiter = " ");


/** Functor for ascending indices sorting i.e., sort the indices based on the corresponding 
 *  values.  
 */
template<typename T>
struct indirect_less_than {
      T *array;
      indirect_less_than(T *value) : array(value) {}
      bool operator() (const int &a, const int &b) {return array[a] < array[b];}
};


/** Functor for descending indices sorting i.e., sort the indices based on the corresponding 
 *  values.  
 */
template<typename T>
struct indirect_greater_than {
      T *array;
      indirect_greater_than(T *value) : array(value) {}
      bool operator() (const int &a, const int &b) {return array[a] > array[b];}
};


/** Write array/vector data into file. Each element occupies a line.
 *
 *  @param fn [read] The file to store the data
 *  @param m [read] Number of elements in data
 *  @param v [read] The data
 */
template <class T>
bool WriteFile(const std::string& fn, const int& m, const T& v)
{
        using namespace std;
        ofstream fp(fn.c_str());
        
        if(not fp.good())
        {
                std::cout << "WriteFile(): Unable to open file (" + fn + ")" << std::endl;
                return false;
        }
        
        for(int i=0; i<m; i++)
                fp << v[i] << std::endl;
        
        fp.close();
        return true;
}


#endif
