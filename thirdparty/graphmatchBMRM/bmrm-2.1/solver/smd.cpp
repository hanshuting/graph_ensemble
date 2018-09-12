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

#ifndef _SMD_CPP_
#define _SMD_CPP_

#include "common.hpp"
#include "timer.hpp"
#include "info.hpp"
#include "configuration.hpp"

#include "smd.hpp"

#include <fstream>
#include <sstream>


/*  Constructor
 
 */
CSMD::CSMD(CModel *model, CLoss *loss)
   : CSolver(model, loss)
{
    Scalar scalar05 = 0.5;
    int numOfHyperplane;
    int dimOfW;
    TheMatrix &w = _model->GetW();
   
    // allocate memory for all attributes
    w.Shape(numOfHyperplane, dimOfW);
    cerr << numOfHyperplane << ";" << dimOfW;
    if(numOfHyperplane!=1) {
        printf("%d",numOfHyperplane);
        throw CBMRMException("number of hyperplanes must be 1!", "SMD::Train()");
    }
    G= new TheMatrix(dimOfW, 1);
    v= new TheMatrix(dimOfW, 1);
    dv= new TheMatrix(dimOfW, 1);
    eta = new TheMatrix(dimOfW, 1);
    temp= new TheMatrix(dimOfW, 1);
    const05 = new TheMatrix(dimOfW, 1);
    
    for(int i=0;i<dimOfW;i++)
        const05->Set(i,scalar05);
    
    // (default) values
    estart = 1.0;
    smdmu = 0.0;
    smdlambda = 0.0;
    tau = 0.0;
    simple_decay = 0;
    project = 0;

    
    accumulate_num  = 0;
    accumulated_obj = 0.0;

}


/** Calculates product of Hessian and v  using finite difference
 *  approximation given parameter w, corresponding gradient g.
 *  The result is stored in hv. temp is a vector used for storing
 *  temporary results which must be preallocated and of same size as g and w.
*/
void CSMD::CalcualeHessianVectorProduct(TheMatrix * w,TheMatrix * g, TheMatrix * v, TheMatrix * temp,TheMatrix *hv) {
   
  
    Scalar epsilon = pow(10.0,-10.0);
    Scalar epsilon1t = pow(10.0,10.0); // 1 / epsilon
    Scalar obj;
    Scalar loss;
    temp->Assign(*v);

    temp->Scale(epsilon); // temp = eps * v 
    w->Add(*temp); //w = w + eps * v
   
    Evaluate(*_loss,loss,obj,*hv);// hv = g(w + eps *v
    w->Minus(*temp);  // reset w
    hv->Minus(*g); // hv = g(w + eps *v) - g(w)
    hv->Scale(epsilon1t);  // hv = (g(w + eps *v) - g(w)) / eps
}


/**  Destructor
 */
CSMD::~CSMD()
{
    delete G;
    delete v;
    delete dv;
    delete eta;
    delete temp;
    delete const05;
     
}



/** Calculates gradient a (with regularized part), loss and objective value
 * given parameter w
 */
void CSMD::Evaluate(CLoss & lossFunction,Scalar & loss,Scalar & obj,TheMatrix & a) {
    loss =0.0;
    lossFunction.ComputeLossAndGradient(loss,a);
    obj =0.0;
    TheMatrix * w = &_model->GetW();
    w->Norm2(obj);
    a.ScaleAdd(lambda,*w);
    obj = obj*obj /2 * lambda;
    
    obj = loss + obj;
}


/**
 * applies max(v_i,max_const) to all elements v_i of vec
*/
void CSMD::maxAllElements(TheMatrix & vec,Scalar max_const) {
    int dim = vec.Length();
    for(int i=0;i<dim;i++) {
        Scalar val;
        vec.Get(i,val);
        if(val<max_const)
            vec.Set(i,max_const);
    }
        

}

/**  Start training/learning a model w.r.t. the loss object (and the data supplied to it).
 *
 */
void CSMD::Train()
{
    
   
    int numOfHyperplane;
    int dimOfW;
 
    //int pos ;
    int    iter = 0;              // iteration count
    Scalar obj = 0.0;            // objective function value        
    Scalar loss = 0.0;            // loss function value        
    
    //double old_obj;                // objective value of last iteration
    CTimer totalTime;             // total runtime of the training
    CTimer lossAndGradientTime;   // time for loss and gradient computation

  
    TheMatrix &w = _model->GetW();
    
  
  
   // check and confirm program parameters from CInfo and/or Configuration file
    ConfirmProgramParameters();


    // retrieve weight vector and gradient dimensionality
    w.Shape(numOfHyperplane, dimOfW);
   


   // start training
   totalTime.Start();
   
   
   v->Zero();
   
    for(int i=0;i<dimOfW;i++)
        eta->Set(i,estart);
    cerr << estart << "eta";
   
   lossAndGradientTime.Start();
   Evaluate(*_loss,loss,obj,*G);
   if (smdmu > 0) { 
       CalcualeHessianVectorProduct(&w, G, v, temp, dv);
   }
   lossAndGradientTime.Stop();
   printf("it %d obj %.10e loss %.10e  time %e \n",iter,obj,loss,totalTime.CPUTotal());
   accumulated_obj=obj;
   while(1) {
       iter++;
       if (smdmu > 0) {
           temp->Assign(*G);
           temp->ElementWiseMult(*v); 
           temp->Scale( - smdmu);  
           temp->ScaleAdd(2.0,*const05);
           
           maxAllElements(*temp,0.5);
           eta->ElementWiseMult(*temp);
      }
       // tau stepsize updata
      if (tau > 0 && iter>1 && simple_decay==0) {
          
          double multiplier =  (tau + iter -2) / (tau + iter -1 );
          eta->Scale(multiplier);
      }
      // 1 /t decay
      if (iter>1 && simple_decay==1) {
          
          double multiplier =  ((double) iter -1) / ((double) iter );
          eta->Scale(multiplier);
      }
      // weight update
      temp->Assign(*G);
      temp->ElementWiseMult(*eta);  
      w.Minus(*temp);
      
      if(project!=0.0) {
          // projection step
          Scalar wlength;
          w.Norm2(wlength); 
          if (wlength>project) {
              printf("scaling from %f to %f\n",wlength,project);
              w.Scale(project / wlength); 
          }
      }
  
      // v    update
      if (smdmu >0) {
          temp->Assign(*G);
          temp->ScaleAdd(smdlambda,*dv); 
          temp->ElementWiseMult(*eta);
          v->Scale(smdlambda); 
          v->Minus(*temp);
      }
      
      /*  gradient */
      lossAndGradientTime.Start();
      Evaluate(*_loss,loss,obj,*G);
      if (smdmu > 0) { 
          CalcualeHessianVectorProduct( &w, G, v, temp, dv);
      }
      lossAndGradientTime.Stop();

      // report loss / objective value / time
      totalTime.Stop();
      double acc_count = 0; // number of accumultion cycles
      double frac_acc_count = 1; // number of minibatches after last accumulation
      if(accumulate_num>0) { // only if we do accumulation
          acc_count = ((double) (iter+1)) / accumulate_num;
          frac_acc_count = (acc_count - ceilMinus(acc_count))*accumulate_num;
      }
     
      // print results of all minibatches of first accumulation cycle
      if(acc_count<1)
          printf("it %d obj %.10e loss %.10e  time %e %e %e\n",iter,obj,loss,totalTime.CPUTotal(),acc_count,frac_acc_count);

      if(frac_acc_count<1.0) { 
          // mini batch is partially in new accumulation cycle
          // split objective value up
          accumulated_obj+= (1.0 - frac_acc_count) * obj;
          printf("ACCUMULATED it %d obj %.10e  time %e \n",(int) ceilMinus(acc_count),accumulated_obj,totalTime.CPUTotal());
          accumulated_obj = frac_acc_count * obj;
      } 
      else {
          accumulated_obj+=obj;
      }
      totalTime.Start();

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


}

 


/**   Confirm/synchronize program parameters set in CInfo and/or Configuration.
 */
void CSMD::ConfirmProgramParameters()
{
   
   Configuration &config = Configuration::GetInstance();  // make sure configuration file is read before this!

   if(config.IsSet("LBFGS.maxNumOfIter")) {
      maxNumOfIter = config.GetInt("SMD.maxNumOfIter");
      if(maxNumOfIter < 0)
          throw CBMRMException("SMD.maxNumOfIter must be > 0\n","SMD::ConfirmProgramParameters()");
   }
   

   if(config.IsSet("SMD.lambda"))           {
      lambda = config.GetDouble("SMD.lambda");
      if(lambda <= 0)
          throw CBMRMException("SMD.lambda must be > 0\n","SMD:ConfirmProgramParameters()");
   }


   if(config.IsSet("SMD.estart"))           {
      estart = config.GetDouble("SMD.estart");
      if(estart <= 0)
          throw CBMRMException("SMD.estart must be > 0\n","SMD::ConfirmProgramParameters()");
   }


   if(config.IsSet("SMD.tau"))           {
       
       tau= config.GetDouble("SMD.tau");
       if(tau < 0)
          throw CBMRMException("SMD.tau must be >= 0\n","SMD::ConfirmProgramParameters()");
   }

   if(config.IsSet("SMD.smdmu"))           {
      smdmu = config.GetDouble("SMD.smdmu");
      if(smdmu < 0)
          throw CBMRMException("SMD.smdmu must be >= 0\n","SMD::ConfirmProgramParameters()");
   }  

   if(config.IsSet("SMD.smdlambda"))           {
      smdlambda = config.GetDouble("SMD.smdlambda");
      if(smdlambda < 0)
          throw CBMRMException("SMD.smdlambda must be >= 0\n","SMD::ConfirmProgramParameters()");
   }

   if(config.IsSet("SMD.simple_decay"))           {
      simple_decay = config.GetInt("SMD.simple_decay");
      if(simple_decay!=0 && simple_decay!=1)
          throw CBMRMException("SMD.simple_decay must be 0 or 1\n","SMD::ConfirmProgramParameters()");
   }


   if(config.IsSet("SMD.project"))           {
      project = config.GetDouble("SMD.project");
      if(project < 0)
          throw CBMRMException("SMD.project must be >= 0\n","SMD::ConfirmProgramParameters()");
   }
   
   if(config.IsSet("Online.minibatches"))           {
       accumulate_num = config.GetDouble("Online.minibatches");
     
   }
   else {
       accumulate_num = 0;
   }

}



#endif

