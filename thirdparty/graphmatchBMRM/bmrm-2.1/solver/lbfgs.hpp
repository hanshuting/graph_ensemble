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

#ifndef _LBFGS_HPP_
#define _LBFGS_HPP_

#include "common.hpp"
#include "solver.hpp"
#include "info.hpp"
#include "loss.hpp"
#include "sml.hpp"   // the matrix library wrapper


/**   Class for LBFGS solver. At the moment a simple
 *    halfing line search is used
 */
class CLBFGS : public CSolver 
{
public:
    // Constructors
    CLBFGS(CModel *model, CLoss *loss);

    // Destructor
    virtual ~CLBFGS();
    
    // Methods
    virtual void Train();
    int INITWITHFIRST;

protected:
    
    static int MAXBACKTRACK;
   
    virtual void ConfirmProgramParameters();
    virtual void ReleaseBuffer();
    virtual void AllocateBuffer(int dim);
    virtual void Evaluate(CLoss &,TheMatrix & w,Scalar & loss,Scalar & obj,TheMatrix & a);
    virtual void BMult(TheMatrix & p);
    virtual void AddToBuffer(TheMatrix & s,TheMatrix & y);
      
    int    maxNumOfIter;      // maximum number of iteration 
     
    Scalar lambda;            // regularization constant
    
    int bufferSize;  
    int actual_bsize; // actual buffer size
    int b_start; // location of first block
    TheMatrix * * S;
    TheMatrix * * Y;
    Scalar * yy;
    Scalar * sy;

};

#endif

