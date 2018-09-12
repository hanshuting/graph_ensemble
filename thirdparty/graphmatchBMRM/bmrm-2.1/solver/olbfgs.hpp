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
 * Created: (73/02/2008) 
 *
 */

#ifndef _OLBFGS_HPP_
#define _OLBFGS_HPP_

#include "common.hpp"
#include "solver.hpp"
#include "info.hpp"
#include "loss.hpp"
#include "lbfgs.hpp"
#include "sml.hpp"   // the matrix library wrapper


/**   Class for oLBFGS solver.
 */
class COLBFGS : public CLBFGS
{
public:
    // Constructors
    COLBFGS(CModel *model, CLoss *loss);

    // Destructor
    virtual ~COLBFGS();
    
    // Methods
    virtual void Train();

protected:
    
    
    double ceilMinus(double v) {
        return ceil(v-1.0);
    }

    virtual void ReleaseBufferSMD();
    virtual void AllocateBufferSMD(int dim);
    //virtual void Evaluate(CLoss & lossFunction,Scalar & loss,Scalar & obj,TheMatrix & a);
    
    void AddToBufferSMD(TheMatrix & s,TheMatrix & ds, TheMatrix & y,TheMatrix & dy);
    void BMultSMD(TheMatrix & p,TheMatrix *dp);
    void CalcualeHessianVectorProduct(TheMatrix * w,TheMatrix * g, TheMatrix * v, TheMatrix * temp,TheMatrix *hv);
    virtual void ConfirmProgramParameters();
   
   
    TheMatrix * v; // SMD v vector
    TheMatrix * dv;  // derivative of v
    TheMatrix * dv2; // derivative of v after update
  

    TheMatrix * * dS; // buffer for derivatives of S
    TheMatrix * * dY; // buffer for derivatives of Y
    
  
    Scalar * dsy; // stores <s,y>
    Scalar * dyy;   // stores <dy,y>
    Scalar * sdy; // stores <s,dy>
  

    // number of batches for which the objective value is accumulated
    // may be non integer for partial batches / if 0 then nothing is
    // accumulated
    Scalar accumulate_num;
    Scalar accumulated_obj; // accumulated ojective 

    //int bufferSize; // buffer size
    //Scalar lambda; // trust region parameter
    Scalar estart; // start learning rate
    Scalar tau; // learning rate decay parameter (if 0 then no decay)

    Scalar smdmu ; // SMD parameter mu (0 -> on SMD used)
    Scalar smdlambda ; // SMD parameter lambda / normal value
    Scalar xlambda ; // SMD parameter lambda / for next iteration
    int usep; // usep1 = 1: use <p,v> instead of <g,v>

    Scalar secant; // target value for trust region adaption
    Scalar tradapt; // speed of changing trust region parameter

    Scalar trpara; // trust-region parameter
  

};

#endif

