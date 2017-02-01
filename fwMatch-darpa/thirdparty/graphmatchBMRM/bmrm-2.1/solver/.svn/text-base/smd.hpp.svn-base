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

#ifndef _SMD_HPP_
#define _SMD_HPP_

#include "common.hpp"
#include "solver.hpp"
#include "info.hpp"
#include "loss.hpp"
#include "sml.hpp"   // the matrix library wrapper


/**   Class for SMD  solver.
 */
class CSMD : public CSolver 
{
public:
    // Constructors
    CSMD(CModel *model, CLoss *loss);

    // Destructor
    virtual ~CSMD();
    
    // Methods
    virtual void Train();

protected:
    
    
    double ceilMinus(double v) {
        return ceil(v-1.0);
    }
   
    virtual void ConfirmProgramParameters();
   
    void Evaluate(CLoss &,Scalar & loss,Scalar & obj,TheMatrix & a);
    void CalcualeHessianVectorProduct(TheMatrix * w,TheMatrix * g, TheMatrix * v, TheMatrix * temp,TheMatrix *hv);
    int    maxNumOfIter;      // maximum number of iterations 
     
    Scalar lambda;            // regularization constant

    /** element-wise maximum operation */
    static void maxAllElements(TheMatrix & vec,Scalar max_const );

    TheMatrix * G; // for storing gradient
    TheMatrix * v;  // SMD v vector
    TheMatrix * dv;  // derivative of v
    TheMatrix * eta;  // learning rates
    TheMatrix * temp;  
    TheMatrix * const05; // constant vector with all entries 0.5

    
    Scalar estart; // start learning rate
    Scalar smdmu ; // SMD parameter mu (if 0 SGD is used)
    Scalar smdlambda ; // SMD parameter lambda
    Scalar tau; // learning rate decay parameter (if 0 then no decay)
    Scalar project; // projection to sphere of length project (0=no projection)
    
    int simple_decay; // if 1 then 1/t deacy is used else tau/(t+tau)


    
    // number of batches for which the objective value is accumulated
    // may be non integer for partial batches / if 0 then nothing is
    // accumulated
    double accumulate_num;
    double accumulated_obj; // accumulated ojective 

};

#endif

