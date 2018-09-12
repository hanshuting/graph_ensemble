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

#ifndef _OLBFGS_CPP_
#define _OLBFGS_CPP_

#include "common.hpp"
#include "olbfgs.hpp"
#include "timer.hpp"
#include "info.hpp"
#include "configuration.hpp"


#include <fstream>
#include <sstream>

// asdf 


/*  Constructor
 
 */
COLBFGS::COLBFGS(CModel *model, CLoss *loss)
   : CLBFGS(model, loss)
{
    
   // set private members (default) values
  
  estart=1;
  tau =1000000;
  lambda=0;
  smdmu =0;
  smdlambda=0;
  secant=0;
  tradapt=0;
  trpara =0.0;
  dS = 0;
  dY = 0;
  dsy = 0;
  sdy = 0;
  usep  = 0;

  
  accumulate_num  = 0;
  accumulated_obj = 0.0;
}


/**  Destructor
 */
COLBFGS::~COLBFGS()
{
   // dont't destruct the loss object! 
   // let the main program does it.

    ReleaseBuffer();
    
}



/** Releases allocated memory for buffer
 *  including SMD attributes
 */
void COLBFGS::ReleaseBufferSMD() {
    ReleaseBuffer(); // also delete LBFGS variables
    delete v;
    delete dv;
    delete dv2;
    
    if(dS) {
        for(int i=0;i<bufferSize;i++)
            delete dS[i];
        delete [] dS;
    } 
    if(dY) {
        for(int i=0;i<bufferSize;i++)
            delete dY[i];
        delete [] dY ;
    }
    if(dsy)
        delete []dsy;
    if(dyy)
        delete []dyy;
    if(sdy)
        delete []sdy;
    dS = 0;
    dY = 0;
    dyy =0;
    dsy = 0;
    sdy = 0;
    actual_bsize = 0;
}

/** Allocates memory for buffer
 *  including SMD attributes
 */
void COLBFGS::AllocateBufferSMD(int dim) {
    ReleaseBufferSMD();
    AllocateBuffer(dim);
    dS = new TheMatrix * [bufferSize];
    dY = new TheMatrix * [bufferSize];
    int row = dim;
    int column = 1;
    for(int i=0;i<bufferSize;i++) {
        dS[i] = new TheMatrix(row, column);
        dY[i] = new TheMatrix(row, column)  ;
    }
    v =  new TheMatrix(row, column);
    dv=  new TheMatrix(row, column);
    dv2 =  new TheMatrix(row, column);
    dsy = new Scalar [bufferSize];
    dyy = new Scalar [bufferSize];
    sdy = new Scalar [bufferSize];
    
}


/** Adds parameter displacement s, its derivative ds and 
    their corresponding change in gradient space,  y and dy, to buffer
*/
void COLBFGS::AddToBufferSMD(TheMatrix & s,TheMatrix & ds, TheMatrix & y,TheMatrix & dy) {
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
    dS[pos]->Assign(ds);
    dY[pos]->Assign(dy);
    S[pos]->Dot(*Y[pos],sy[pos]);
    Y[pos]->Dot(*Y[pos],yy[pos]);
    dY[pos]->Dot(*Y[pos],dyy[pos]);
    dS[pos]->Dot(*Y[pos],dsy[pos]);
    S[pos]->Dot(*dY[pos],sdy[pos]);
    
}


 
/** Calculates product of Hessian and v  using finite difference
 *  approximation given parameter w, corresponding gradient g.
 *  The result is stored in hv. temp is a vector used for storing
 *  temporary results which must be preallocated and of same size as g and w.
*/
void COLBFGS::CalcualeHessianVectorProduct(TheMatrix * w,TheMatrix * g, TheMatrix * v, TheMatrix * temp,TheMatrix *hv) {
    Scalar epsilon = pow(10.0,-10.0);
    Scalar epsilon1t = pow(10.0,10.0);// 1 / epsilon
    Scalar obj;
    Scalar loss;
    temp->Assign(*v);

    temp->Scale(epsilon); // temp = eps * v 
    w->Add(*temp); //w = w + eps * v
    Evaluate(*_loss,*w,loss,obj,*hv); //hv = g(w + eps *v) 
    w->Minus(*temp); // reset w
    hv->Minus(*g); // hv = g(w + eps *v) - g(w)
    hv->Scale(epsilon1t);  // hv = (g(w + eps *v) - g(w)) / eps
}




/** Calulates implizit product of B
    (approximation of inverse Hessian) wit p and dp.
    In case of dp==0 only Bp is calculated
 */
void COLBFGS::BMultSMD(TheMatrix & p,TheMatrix  *dp) { 
    Scalar db0;
    Scalar b0;

    int pos;
    Scalar * alpha = new Scalar [bufferSize];
    Scalar * dalpha = new Scalar [bufferSize];

    //Scalar beta;

   
    for (int i=actual_bsize-1;i>=0;i--) {
        pos = (b_start + i) % bufferSize;
        Scalar powlambda = pow(xlambda,actual_bsize-i);
        Scalar sp;
        Scalar rho = 1.0/sy[pos];
        S[pos]->Dot(p,sp);
        alpha[i] = rho *sp;
        
        
        if(dp!=0) {
            Scalar dsp =0;
            Scalar sdp =0;
            dS[pos]->Dot(p,dsp);
            S[pos]->Dot(*dp,sdp);
            dalpha[i]= -  rho * rho * powlambda * sp* (sdy[pos] + dsy[pos]) + rho* (powlambda * dsp + sdp);
            
        }
        p.ScaleAdd(-alpha[i],*Y[pos]);
        if (dp!=0) {
            dp->ScaleAdd(- dalpha[i],*Y[pos]);
            dp->ScaleAdd( -alpha[i]* powlambda,*dY[pos]);
        }
    }
      
    // averaging b0 and db0 values over whole buffer 
    b0 = 0.0;
    db0 = 0.0;
    for (int i=0;i<actual_bsize;i++) {
        Scalar powlambda = pow(xlambda,actual_bsize-i);
        pos = (b_start + i) % bufferSize;
        b0+=  sy[pos]/ yy[pos];
        if (dp!=0) {
            db0+=  powlambda * (-2 * sy[pos] *dyy[pos] / (yy[pos]*yy[pos]) + (dsy[pos] + sdy[pos]) / yy[pos]);
        }
       
      }
    b0= b0 / ((Scalar) actual_bsize);
    db0= db0 /  ((Scalar) actual_bsize);
    
    if (dp!=0) {
       dp->Scale(b0); 
       dp->ScaleAdd(db0,p); 
    }
    p.Scale(b0);
    for (int i=0;i<actual_bsize;i++) {
        Scalar dbeta = 0.0;
        Scalar beta;
        Scalar powlambda = pow(xlambda,actual_bsize-i);
        pos = (b_start + i) %  bufferSize;
        Scalar rho = 1.0/sy[pos];
        Scalar yp;
        Y[pos]->Dot(p,yp);
        beta = rho * yp;
        if(dp!=0) {
            Scalar ydp;
            Scalar dyp;
            dY[pos]->Dot(p,dyp);
            Y[pos]->Dot(*dp,ydp);
            dbeta = - powlambda * rho * rho * yp * (sdy[pos] +dsy[pos]) +  rho *(dyp * powlambda  + ydp); 
        }
        
        p.ScaleAdd(alpha[i]-beta,*S [pos]);
        if (dp!=0 ) {
            dp->ScaleAdd(dalpha[i]-dbeta,*S[pos]); 
            dp->ScaleAdd((alpha[i]-beta)* powlambda,*dS[pos]); 
            
      }

    }
    delete [] dalpha;
    delete [] alpha;
}



/**  Start training/learning a model w.r.t. the loss object (and the data supplied to it).
 *
 */
void COLBFGS::Train()
{
  
   
    int numOfHyperplane;
    int dimOfW;

    b_start =0;
    actual_bsize=0;
    //int pos ;
    int    iter = 0;              // iteration count
    //int    evals = 0;              // evaluation counts
    Scalar obj = 0.0;            // objective function value        
    Scalar loss = 0.0;            // loss function value        
    
    //Scalar old_obj;                // objective value of last iteration
    CTimer totalTime;             // total runtime of the training
    CTimer lossAndGradientTime;   // time for loss and gradient computation

  
    TheMatrix &w = _model->GetW();
    
  

   // check and confirm program parameters from CInfo and/or Configuration file
   ConfirmProgramParameters();

   // retrieve weight vector and gradient dimensionality
   w.Shape(numOfHyperplane, dimOfW);
   if(numOfHyperplane!=1) {
       printf("%d",numOfHyperplane);
       throw CBMRMException("number of hyperplanes must be 1!", "COLBFGS::Train()");
   }

   
   AllocateBufferSMD(dimOfW);
   TheMatrix * p =new TheMatrix(dimOfW, 1); // search direction
   TheMatrix* dp = new TheMatrix(dimOfW, 1); // derivative of search direction
   
   TheMatrix * G = new TheMatrix(dimOfW, 1); // gradient
   TheMatrix * g2 = new TheMatrix(dimOfW, 1); // gradient after update
   TheMatrix * temp = new TheMatrix(dimOfW, 1); 

   Scalar       first_stepsize = 0.0000001; // used for the gradient descent step
   Scalar e = estart;
   xlambda = smdlambda;
   // start training
   totalTime.Start();
   
   lossAndGradientTime.Start();
   Evaluate(*_loss,w,loss,obj,*G);
   if (smdmu > 0) { 
       CalcualeHessianVectorProduct(&w, G, v, temp, dv);
   }
   lossAndGradientTime.Stop();
   //printf("it %d obj %.10e loss %.10e  time %e \n",iter,obj,loss,totalTime.CPUTotal());
   totalTime.Start();
   accumulated_obj=obj; 

   v->Zero();
   p->Assign(*G);;
   p->Scale(-first_stepsize); 
   w.Add(*p);
   if (smdmu > 0) {
       dp->Assign(*p);

       v->Add(*dp);
   }

   /* Second gradient */
   lossAndGradientTime.Start();
   Evaluate(*_loss,w,loss,obj,*g2);
   if (smdmu > 0) {
       CalcualeHessianVectorProduct(&w, g2, v, temp, dv2);
       dv2->ScaleAdd(-xlambda,*dv);
       dv2->ScaleAdd(trpara,*dp);      
   }
   lossAndGradientTime.Stop();
   g2->Minus(*G);
   g2->ScaleAdd(trpara,*p);
   
   if(smdmu >0)
       AddToBufferSMD(*p,*dp,*g2,*dv2);
   else
       AddToBuffer(*p,*g2);
   
   while(1) {
       iter++;
       lossAndGradientTime.Start();
      
       Evaluate(*_loss,w,loss,obj,*G);
       
       if (smdmu > 0) { 
           CalcualeHessianVectorProduct(&w, G, v, temp, dv);
       }
       lossAndGradientTime.Stop();
          
       // report loss / objective value / time
       totalTime.Stop();
       Scalar acc_count = 0; // number of accumultion cycles
       Scalar frac_acc_count = 1; // number of minibatches after last accumulation
 
      if(accumulate_num>0) { // only if we do accumulation
           acc_count = ((Scalar) (iter+1)) / accumulate_num;
           frac_acc_count = (acc_count - ceilMinus(acc_count))*accumulate_num;
       
       }
       // print results of all minibatches of first accumulation cycle
       if(acc_count<1.0)
           printf("it %d obj %.10e loss %.10e  time %e %e %e\n",iter,obj,loss,totalTime.CPUTotal(),acc_count,frac_acc_count);
       
       if(frac_acc_count<=1.0) {
           // mini batch is partially in new accumulation cycle
           // split objective value upnext accumulation
           accumulated_obj+= (1.0 - frac_acc_count) * obj;
           printf("ACCUMULATED it %d obj %.10e  time %e %e %e \n",(int) ceilMinus(acc_count),accumulated_obj,totalTime.CPUTotal(),acc_count,frac_acc_count);
           accumulated_obj = frac_acc_count * obj;
       } 
       else {
           accumulated_obj+=obj;
       }
       totalTime.Start();

       p->Assign(*G);
       if (smdmu >0 ) {
           dp->Assign(*dv);
           dp->Scale(xlambda); 
          
       }
      
       BMultSMD(*p,dp);
       

       if (smdmu>0) {

           Scalar stepchange = 0.0;
            // p instead of gg gradient
           if (usep) {
               p->Dot(*v,stepchange);
           }
           else {
              G->Dot(*v,stepchange);

           }

           
           stepchange = stepchange *smdmu;
           xlambda = 0.0;
           if (stepchange > 0.5) {
               stepchange=0.5;
           }
           else {
               if (stepchange <-1)
                   stepchange = -1;
               else
                   xlambda = smdlambda;
           }
           e = e * (1 - stepchange);
           /* if (iter % 1000 ==0)
              printf("e %f\n",e);*/
       }
       // gain decay update
      if (tau > 0 && iter>1) {
          Scalar multiplier =  (tau + iter -2) / (tau + iter -1 );
          e*=multiplier;
      }
     
      p->Scale(-e);
      dp->Scale(-e);
      
      // weight update
      w.Add(*p);
      if (smdmu > 0) {
          dp->Add(*p);
          v->Scale(xlambda);
          v->Add(*dp);
      }
      /* Second gradient */
      lossAndGradientTime.Start();
     
      Evaluate(*_loss,w,loss,obj,*g2);
      
      
      if (smdmu > 0) {
          CalcualeHessianVectorProduct(&w, g2, v, temp, dv2);
          dv2->ScaleAdd(-xlambda,*dv);
          dv2->ScaleAdd(trpara,*dp);      
      }
      lossAndGradientTime.Stop();
      
      g2->ScaleAdd(-1.0,*G);
      g2->ScaleAdd(trpara,*p);
      if(tradapt > 0) {
          
          temp->Assign(*g2); // temp is current y
          
          BMultSMD(*temp,0);
          Scalar norm_p;
          Scalar norm_temp;
          p->Norm2(norm_p); 
          temp->Minus(*p); // temp = By -s
          temp->Norm2(norm_temp); 
          Scalar multiplier = 1 + tradapt * (norm_temp/norm_p - secant);
          if (multiplier > 2)
              multiplier = 2;
          if(multiplier < 0.5)
              multiplier = 0.5;           
          trpara*=multiplier;
          /*if(iter % 1000 ==0)
            printf("trust r. %f; %f\n",norm_temp/norm_p,trpara); */
     
      }
      if(smdmu >0)
          AddToBufferSMD(*p,*dp,*g2,*dv2);
      else
          AddToBuffer(*p,*g2);
   
    
      if(iter >= maxNumOfIter ) 
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

   delete temp;
   delete p;
   delete dp;
   delete G;
   delete g2;
   

}

 /**   Confirm/synchronize program parameters set in CInfo and/or Configuration.
 */
void COLBFGS::ConfirmProgramParameters()
{
   
   Configuration &config = Configuration::GetInstance();  // make sure configuration file is read before this!
   if(config.IsSet("OLBFGS.maxNumOfIter")) {
      maxNumOfIter = config.GetInt("OLBFGS.maxNumOfIter");
      
      if(maxNumOfIter < 0)
          throw CBMRMException("OLBFGS.maxNumOfIter must be > 0\n","COLBFGS::ConfirmProgramParameters()");
   }
   
   if(config.IsSet("OLBFGS.buffersize")) {
      bufferSize= config.GetInt("OLBFGS.buffersize");  
      if(bufferSize <= 0)
          throw CBMRMException("OLBFGS.buffersize must be > 0\n","COLBFGS::ConfirmProgramParameters()");
   }



   if(config.IsSet("OLBFGS.lambda"))           {
      lambda = config.GetDouble("OLBFGS.lambda");
      if(lambda <= 0)
          throw CBMRMException("OLBFGS.lambda must be > 0\n","OLBFGS:ConfirmProgramParameters()");
   }


   if(config.IsSet("OLBFGS.estart"))           {
      estart = config.GetDouble("OLBFGS.estart");
      if(estart <= 0)
          throw CBMRMException("OLBFGS.estart must be > 0\n","OLBFGS::ConfirmProgramParameters()");
   }


   if(config.IsSet("OLBFGS.tau"))           {
       
       tau= config.GetDouble("OLBFGS.tau");
       if(tau < 0)
          throw CBMRMException("OLBFGS.tau must be >= 0\n","OLBFGS::ConfirmProgramParameters()");
   }
   if(config.IsSet("OLBFGS.trpara"))           {
       
       trpara= config.GetDouble("OLBFGS.trpara");
       if(trpara < 0)
          throw CBMRMException("OLBFGS.trpara must be >= 0\n","OLBFGS::ConfirmProgramParameters()");
   }

   if(config.IsSet("OLBFGS.smdmu"))           {
      smdmu = config.GetDouble("OLBFGS.smdmu");
      if(smdmu < 0)
          throw CBMRMException("OLBFGS.smdmu must be >= 0\n","OLBFGS::ConfirmProgramParameters()");
   }  

   if(config.IsSet("OLBFGS.smdlambda"))           {
      smdlambda = config.GetDouble("OLBFGS.smdlambda");
      if(smdlambda < 0)
          throw CBMRMException("OLBFGS.smdlambda must be >= 0\n","OLBFGS::ConfirmProgramParameters()");
   }


   if(config.IsSet("OLBFGS.usep"))           {
      usep = config.GetInt("OLBFGS.usep");
      if(usep!=0 and usep!=1)
          throw CBMRMException("OLBFGS.usep must be 0 or 1\n","OLBFGS::ConfirmProgramParameters()");
   }
   
   if(config.IsSet("OLBFGS.secant"))           {
      secant = config.GetDouble("OLBFGS.secant");
   }  
   if(config.IsSet("OLBFGS.tradapt"))           {
      tradapt = config.GetDouble("OLBFGS.tradapt");
      if(tradapt < 0)
          throw CBMRMException("OLBFGS.tradapt must be >= 0\n","OLBFGS::ConfirmProgramParameters()");
   }
    
   if(config.IsSet("Online.minibatches"))           {
       accumulate_num = config.GetDouble("Online.minibatches");
     
   }
   else {
       accumulate_num = 0;
   }

}



#endif

