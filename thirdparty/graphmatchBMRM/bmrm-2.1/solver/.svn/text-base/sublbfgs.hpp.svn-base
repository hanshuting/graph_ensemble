/* Copyright (c) 2007, National ICT Australia 
 * All rights reserved. 
 * 
 * The contents of this file are subject to the Mozilla Public License 
 * Version 1.1 (the "License"); you may not use this file except in 
 * compliance with the License. You may obtain a copy of the License at 
 * http://www.mozilla.org/MPL/ 
 * 
 * Software distributed under the License is distributed on an "AS IS" 
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the 
 * License for the specific language governing rights and limitations 
 * under the License. 
 * 
 * Authors: Simon Guenter (simon.guenter@nicta.com.au)
 *
 * Created: (03/02/2008) 
 *
 */

#ifndef _SUBLBFGS_HPP_
#define _SUBLBFGS_HPP_

#include "lbfgs.hpp" 

#include "common.hpp"
#include "solver.hpp"
//#include "bmrminnersolver.hpp"
#include "model.hpp"
#include "loss.hpp"

#include "vecdata.hpp"

#include "timer.hpp"
#include <algorithm>
#include <vector>

/** for sorting indices */
class ScalarWithIndex {
public:
    Scalar value;
    int index;
    ScalarWithIndex(Scalar val,int ind): value(val),index(ind) {};
    ScalarWithIndex(): value(0.0),index(-1) {};
    ~ScalarWithIndex() {};
    inline int operator<(const ScalarWithIndex & other) const
    {
        return(value < other.value);
    }
    

    inline void swap(ScalarWithIndex & other) {
        Scalar myvalue = value;
        int myindex = index;
        value = other.value;
        index = other.index;
        other.value = myvalue;
        other.index = myindex;
    }
};

/**   Class for SUBLBFGS solver.
 */
class CSUBLBFGS : public CLBFGS 
{
public:
   
    CSUBLBFGS(CModel * model, CVecData * data);

    // Destructor
    virtual ~CSUBLBFGS();
    
    
    virtual void Train();
    static void sortIndices(std::vector<ScalarWithIndex> & v);
    Scalar lineSearch2(TheMatrix & p,TheMatrix & w,TheMatrix &f,TheMatrix * temp,TheMatrix * beta,TheMatrix * df,vector<ScalarWithIndex> * alpha,vector<int> &hinges);
    Scalar lineSearch(TheMatrix & p,TheMatrix & w,TheMatrix &f,TheMatrix * temp,TheMatrix * beta,TheMatrix * df,vector<ScalarWithIndex> * alpha,vector<int> &hinges);
  
    int getDescentDirection(TheMatrix & p,TheMatrix &g,TheMatrix & w,vector<int> &hinges,vector<double> &hingebetas,TheMatrix * a2, TheMatrix * gnew,TheMatrix * temp1,TheMatrix * temp2, TheMatrix * ga);
    void evaluateFunction(TheMatrix &w, Scalar & loss,Scalar & obj,TheMatrix &g,vector<int> &hinges,vector<double> &hingebetas,TheMatrix & beta,TheMatrix & f,bool calc_f);
    
    void setC(const double& cVal){C = cVal;}

    int USENEWGRADIENT;
    int MAXDESCENTTRY;
    int DEBUGLEVEL;
protected: 
    CVecData *_data;
    CTimer totalTime;             // total runtime of the training
    CTimer lossAndGradientTime;   // time for loss and gradient computation
    CTimer lineSearchTime;   // time for loss and gradient computation
    CTimer descentDirectionTime;   // time for loss and gradient computation
    double C;
    int numElements;
   

};

#endif
