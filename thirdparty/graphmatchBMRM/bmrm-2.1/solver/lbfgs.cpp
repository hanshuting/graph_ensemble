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

#ifndef _LBFGS_CPP_
#define _LBFGS_CPP_

#include "common.hpp"
#include "lbfgs.hpp"
#include "timer.hpp"
#include "info.hpp"
#include "configuration.hpp"


#include <fstream>
#include <sstream>

int CLBFGS::MAXBACKTRACK = 30;


/*  Constructor
 
 */
CLBFGS::CLBFGS(CModel *model, CLoss *loss)
   : CSolver(model, loss)
{
    
   // set private members (default) values
   maxNumOfIter     = 1000;
   lambda           = 1.0;
   bufferSize = 5;
   INITWITHFIRST = 0;

   S = 0;
   Y = 0;
   yy = 0;
   sy = 0;
}


/**  Destructor
 */
CLBFGS::~CLBFGS()
{
   // dont't destruct the loss object! 
   // let the main program does it.

    ReleaseBuffer();
    
}



/** Releases allocated memory for buffer */
void CLBFGS::ReleaseBuffer() {
    if(S) {
        for(int i=0;i<bufferSize;i++)
            delete S[i];
        delete [] S;
    } 
    if(Y) {
        for(int i=0;i<bufferSize;i++)
            delete Y[i];
        delete []Y ;
    }
    if(yy)
        delete []yy;
    if(sy)
        delete []sy;
    S = 0;
    Y = 0;
    yy =0;
    sy = 0;
    actual_bsize = 0;
}

/** Allocates memory for buffer */
void CLBFGS::AllocateBuffer(int dim) {
    ReleaseBuffer();
    S = new TheMatrix * [bufferSize];
    Y = new TheMatrix * [bufferSize];
    int row = dim;
    int column = 1;
    for(int i=0;i<bufferSize;i++) {
        S[i] = new TheMatrix(row, column);
        Y[i] = new TheMatrix(row, column)  ;
    }
    yy = new Scalar [bufferSize];
    sy = new Scalar [bufferSize];
}


/** Adds pair of parameter displacement s and
    corresponding change in gradient y to buffer
*/
void CLBFGS::AddToBuffer(TheMatrix & s,TheMatrix & y) {
    int pos;
    actual_bsize++;
   
    if (actual_bsize>bufferSize) {
        // buffer full, delete last element and add new one
        pos = b_start    ;
        b_start = (b_start + 1)  %bufferSize;
        actual_bsize--;
    }
    else {
        pos = actual_bsize - 1;
    }
    
    S[pos]->Assign(s);
    Y[pos]->Assign(y);
    S[pos]->Dot(*Y[pos],sy[pos]);
    Y[pos]->Dot(*Y[pos],yy[pos]);
    
}


/** Calculates gradient a (with regularized part), loss and objective value
 * given parameter w
 */
void CLBFGS::Evaluate(CLoss & lossFunction,TheMatrix & w,Scalar & loss,Scalar & obj,TheMatrix & a) {
    loss =0.0;
    lossFunction.ComputeLossAndGradient(loss,a);
    obj =0.0;
    w.Norm2(obj);
    a.ScaleAdd(lambda,w);
    obj = obj*obj /2 * lambda;
    obj = loss + obj;

}



/** Calulates implizit product of B
    (approximation of inverse Hessian) and p */
void CLBFGS::BMult(TheMatrix & p) {
    int pos;
    Scalar * alpha = new Scalar [bufferSize];

    Scalar beta;
    for (int i=actual_bsize-1;i>=0;i--) {
        pos = (b_start + i) % bufferSize;
        //pos = b_start;
        S[pos]->Dot(p,alpha[i]);
        alpha[i] = 1.0/ sy[pos] *alpha[i];
        p.ScaleAdd(-alpha[i],*Y[pos]);
    }
      
    
      if(actual_bsize>0) {
          pos = (b_start + actual_bsize -1) % bufferSize;
          if (INITWITHFIRST)
              pos = b_start;
          p.Scale(sy[pos] / yy[pos]);
      }
      for (int i=0;i<actual_bsize;i++) {
          pos = (b_start + i) %  bufferSize;
          
          Y[pos]->Dot(p,beta);
          beta = 1.0 / sy[pos] * beta;
        
          p.ScaleAdd(alpha[i]-beta,*S [pos]);
      }
      delete [] alpha;
}


/**  Start training/learning a model w.r.t. the loss object (and the data supplied to it).
 *
 */
void CLBFGS::Train()
{
    
   
    int numOfHyperplane;
    int dimOfW;

    b_start =0;
    actual_bsize=0;
    //int pos ;
    int    iter = 0;              // iteration count
    int    evals = 0;              // evaluation counts
    Scalar obj = 0.0;            // objective function value        
    Scalar loss = 0.0;            // loss function value        
    
    Scalar old_obj;                // objective value of last iteration
    CTimer totalTime;             // total runtime of the training
    CTimer lossAndGradientTime;   // time for loss and gradient computation

  
    TheMatrix &w = _model->GetW();
    
  
    //CLoss &lossFunction = _model->GetLoss();
   TheMatrix* a = 0;                // gradient
   TheMatrix * p =0;
   TheMatrix* a2 = 0;
  
   // check and confirm program parameters from CInfo and/or Configuration file
   ConfirmProgramParameters();


   // retrieve weight vector and gradient dimensionality
   w.Shape(numOfHyperplane, dimOfW);
   if(numOfHyperplane!=1) {
       printf("%d",numOfHyperplane);
       throw CBMRMException("number of hyperplanes must be 1!", "CLBFGS::Train()");
   }

   AllocateBuffer(dimOfW);
   
   a = new TheMatrix(dimOfW, 1);
   a2 = new TheMatrix(dimOfW, 1);
   p = new TheMatrix(dimOfW,1);
   
   // start training
   totalTime.Start();
   
   // gradient calculation
   lossAndGradientTime.Start();
   Evaluate(*_loss,w,loss,obj,*a);
   lossAndGradientTime.Stop();
       
   while(1) {
      
       iter++;
       evals++;
      
    
      // search direction calculation
      p->Assign(*a);
      BMult(*p);
      
      // line Search (here simply halfing until finished)
      Scalar eta = 1.0;
      old_obj = obj; 
      w.ScaleAdd(-1.0,*p);
      lossAndGradientTime.Start();
      Evaluate(*_loss,w,loss,obj,*a2);
      lossAndGradientTime.Stop();
      int backtrackcnt = 0;
      while(old_obj<=obj) {
          eta = eta * 0.5;
          w.ScaleAdd(eta,*p);
          lossAndGradientTime.Start();
          Evaluate(*_loss,w,loss,obj,*a2);
          lossAndGradientTime.Stop();
          totalTime.Stop();
         
          totalTime.Start();
          backtrackcnt++;
          evals++;
          if(backtrackcnt>MAXBACKTRACK) {
              break;
          }
      }
      if(backtrackcnt>MAXBACKTRACK) {
          printf("max number of backtracks in line search exceeded");
          w.ScaleAdd(eta,*p);
          
          break;
      }
      p->Scale(-eta);
      a->Scale(-1.0);
      a->Add(*a2);

      totalTime.Stop();

      printf("it %d eval %d obj %.10e loss %.10e  time %e \n",iter,evals,obj,loss,totalTime.CPUTotal());
      totalTime.Start();


      /* buffer management */
      AddToBuffer(*p,*a);
      
      a->Assign(*a2);
      
      if(iter >= maxNumOfIter)
      { 
         printf("\nWARNING: program exceeded maximum number of iterations (%d) !\n", maxNumOfIter);
         break;
      }
   }

 
   
   totalTime.Stop();

   // display timing profile
   printf("\nCPU seconds in:\n");
   printf("1. loss and gradient: %8.2f\n", lossAndGradientTime.CPUTotal());
  
   printf("               Total: %8.2f\n", totalTime.CPUTotal());
   printf("\nLegends: \nit: iteration number \neval: number of function evaluations\nobj: objective value\nloss: loss function value\ntime: CPU time used\n\n"); 

   delete a;               
   delete p;
   delete a2;

}

 


/**   Confirm/synchronize program parameters set in CInfo and/or Configuration.
 */
void CLBFGS::ConfirmProgramParameters()
{
   
   Configuration &config = Configuration::GetInstance();  // make sure configuration file is read before this!

   if(config.IsSet("LBFGS.maxNumOfIter")) {
      maxNumOfIter = config.GetInt("LBFGS.maxNumOfIter");
      if(maxNumOfIter < 0)
          throw CBMRMException("LBFGS.maxNumOfIter must be > 0\n","CLBFGS::ConfirmProgramParameters()");
   }
   
   if(config.IsSet("LBFGS.bufferSize")) {
      bufferSize= config.GetInt("LBFGS.bufferSize");  
      if(bufferSize <= 0)
          throw CBMRMException("LBFGS.bufferSize must be > 0\n","CLBFGS::ConfirmProgramParameters()");
   }

   if(config.IsSet("LBFGS.lambda"))           {
      lambda = config.GetDouble("LBFGS.lambda");
      if(lambda <= 0)
	 throw CBMRMException("LBFGS.lambda must be > 0\n","CLBFGS::ConfirmProgramParameters()");
   }



}



#endif

