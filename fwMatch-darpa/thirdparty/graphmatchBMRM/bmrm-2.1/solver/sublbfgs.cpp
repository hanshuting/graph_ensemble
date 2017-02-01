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

#ifndef _SUBLBFGS_CPP_
#define _SUBLBFGS_CPP_

#include "common.hpp"
#include "lbfgs.hpp"
#include "sublbfgs.hpp"

#include "timer.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
//#include "bmrminnersolver.hpp"
#include "loss.hpp"
//#include "bmrminnersolverfactory.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>



/** evaluates loss, objective value (obj) and gradient (g) for
    parameter w using hinge loss. Also returns the element numbers
    for which w is on the hinge (hinges) and the corresponding
    randomly chosen value which was used for the
    gradient g (beta). f is X*w, X being the data matrix.
    If calculate_f==1 then f is recalculated.
*/
void CSUBLBFGS::evaluateFunction(TheMatrix &w, Scalar & loss,Scalar & obj,TheMatrix &g,vector<int> &hinges,vector<double> &hingebetas,TheMatrix & beta,TheMatrix & f,bool calculate_f=1) {

    if(calculate_f) {
        _data->ComputeF(w,f);
	f.ElementWiseMult(_data->labels());
    }
    
    hingebetas.clear();
    hinges.clear();
    
    loss = 0;  
    for (int i=0;i<numElements;i++) {
        Scalar val;
        f.Get(i,val);
        beta.Set(i,0);
        Scalar error = 1.0-val;
        if(error==0 ) {
            hinges.push_back(i);
	    Scalar randomval = rand()/(double(RAND_MAX)+1);
            beta.Set(i, -C*randomval); 
            hingebetas.push_back(randomval);  
        }
        if(error>0) {
            loss+= error;
            beta.Set(i,-C); 
        }
    }
   
   
    loss*= C;
    beta.ElementWiseMult(_data->labels());

    _data->Grad(beta,g);
    g.ScaleAdd(lambda,w); // add regularizer
    w.Norm2(obj);
    obj= 0.5 * lambda * obj*obj; // regularizer
    obj+=loss;
    if(DEBUGLEVEL)
        cout << "number of hinges: " << hinges.size() << endl;
    
}

void CSUBLBFGS::sortIndices(vector<ScalarWithIndex> & v) {
    
    sort(v.begin(), v.end());
    
}

       
      

/*  Constructor
 
 */
CSUBLBFGS::CSUBLBFGS(CModel *model,CVecData * data)
  : CLBFGS(model, 0)
{
   
    _data =  data;
    MAXDESCENTTRY = 1000;
    USENEWGRADIENT = 1;
    DEBUGLEVEL = 0;
    numElements = _data->size();
    C = 1.0/(double)numElements;
}


/**  Destructor
 */
CSUBLBFGS::~CSUBLBFGS() 
{  
}


/**
   returns optimal step size starting from w in the
   direction of p. f=Xw with X being the data matrix.
   beta,df and alpha are temp variables of size of the
   data set which must be initialized by the caller
   This method only does partial sorting starting with
   stepsizes smaller than 2 and doubles this stepsize
   if minimum is not found.
*/
Scalar CSUBLBFGS::lineSearch2(TheMatrix & p,TheMatrix & w,TheMatrix &f,TheMatrix * temp,TheMatrix * beta,TheMatrix * df,vector<ScalarWithIndex> * alpha,vector<int> &hinges)
{

    // 
    hinges.clear();
    //
    
    // start searching with interval that starts before
    // min_step and ends after min_step
    Scalar min_stepsize = -pow(10.0,-10.0); 

     //Time0.Start();
    //TheMatrix * p2 = new TheMatrix(dimOfW,1);
    //p2->Assign(p);
    //p2->Scale(-1.0);
    
    lineSearchTime.Stop();
    lossAndGradientTime.Start();
    _data->ComputeF(p,*df);
    df->ElementWiseMult(_data->labels());
    lossAndGradientTime.Stop();
    lineSearchTime.Start();
    
    Scalar b;
    Scalar h;
    //Time0.Stop();
    //Time1.Start();
    for (int i=0;i<numElements;i++) {
        Scalar val1;
        Scalar val2;
        f.Get(i,val1);
        df->Get(i,val2);
        (*alpha)[i].index = i;
        (*alpha)[i].value = (1- val1 ) / val2;
       
    }
    p.Dot(w,b);
    b*=lambda;
    
    p.Norm2(h);
    h=h*h;
    h*= lambda;
    //Time1.Stop();
    
    
    
    // PRESORTING
    //put all elements smaller than min_stepsize before
    //the elements larger than min_stepsize
    //also stores largest number smaller than min_stepsize
    Scalar maxneg = 0.0; //   largest number smaller than min_stepsiz
    int maxnegindex = -1; // index of largest number smaller than min_stepsiz
    int firstposindex = -1; // all indices smaller than firstposindex are numbers smaller than min_stepsize
    
    // until we find a number larger than min_stepsize
    for(int i=0;i<numElements;i++) {
        Scalar thisval = (*alpha)[i].value;
        //int thisindex = (*alpha)[i].index;
        //store largest number
        if(thisval<min_stepsize) {
            if(thisval>maxneg || maxnegindex ==-1) {
                maxnegindex = i;
                maxneg = thisval;
            }
        }
        else {
            // first number found -> break
            firstposindex = i;
            break;

        }
    }
    // swapping will gives us the desired ordering
    for(int i=firstposindex+1;i<numElements;i++) {
        Scalar thisval = (*alpha)[i].value;
        //int thisindex = (*alpha)[i].index;
        if(thisval<min_stepsize) {
            if(thisval>maxneg || maxnegindex ==-1) {
                maxnegindex = i;
                maxneg = thisval;
            }
            
            (*alpha)[i].swap((*alpha)[firstposindex]);
            
            if(maxneg==thisval) 
                maxnegindex = firstposindex;
            firstposindex++;
                
        }
    }

    //cout << "maxnegindex" << maxnegindex << maxneg;
    if(maxnegindex!=-1) {
        // "move" largest number smaller then min_stepsize
        // one postion before all numbers larger then min_stepsize
        (*alpha)[maxnegindex].swap((*alpha)[firstposindex-1]);
        firstposindex--;
    }
    
    // all numbers before firstposindex are irrelevant and
    // can have any kind of ordering

    if(firstposindex <0)
        firstposindex = 0;
    

    // only sort numbers between min_stepsize and max_stepsize
    // if we don't find optimal step size increase max_stepsize
    // and check again (excluding the already checked alphas)
    Scalar max_stepsize = 2.0; // first guess of upper bound of eta
    int lastposindex = numElements-1;   


    Scalar eta;
    int firstval = 1;
    int lastval =1;
    int first = 1;


    
    while(1) {
        for(;lastposindex>=firstposindex;lastposindex--) {
            Scalar thisval = (*alpha)[lastposindex].value;
            if(thisval<max_stepsize)
                break;
        }
        // as long as there is no element, continue increasing max step size
        while(lastposindex==firstposindex) {
            max_stepsize*=2;
            if(DEBUGLEVEL) 
                cout << "increasing max stepsize to " << max_stepsize << ".\n";
            lastposindex =  numElements-1;
            for(;lastposindex>=firstposindex;lastposindex--) {
                Scalar thisval = (*alpha)[lastposindex].value;
                if(thisval<max_stepsize)
                    break;
            }
        }
        for(int i=lastposindex;i>=firstposindex;i--) {
            
            Scalar thisval = (*alpha)[i].value;
           
            if(thisval>=max_stepsize) {
                (*alpha)[i].swap((*alpha)[lastposindex]);
                lastposindex--;
                
            }
        }
        if(DEBUGLEVEL)
            cout << " first index: " << (firstposindex+1) << " last index:  " <<  lastposindex << "\n";
        
        sort(alpha->begin() + firstposindex,alpha->begin() + lastposindex+1);
  

        int i = firstposindex + 1;
        if(first) {
            first = 0;
            eta = ((*alpha)[i-1].value +  (*alpha)[i].value) / 2.0;
            temp->Assign(f);
            temp->ScaleAdd(eta,*df);
            for (int j=0;j<numElements;j++) {
                Scalar val;
                temp->Get(j,val);
                val = (val<1);
                beta->Set(j,val);
            }
            beta->Dot(*df,eta);
            eta = (C * eta - b) / h; 
        }
        
        
        
        while(i<=lastposindex) {
     
            if(eta>=(*alpha)[i-1].value) {
                if(eta<=(*alpha)[i].value)
                    break;
                Scalar flip = (*alpha)[i].value;
                firstval = i;
                while((*alpha)[i].value==flip) { // for all hinges with same alpha
                    
                    int index = (*alpha)[i].index;
                    Scalar val;
                    Scalar dfval;
                    beta->Get(index,val);
                    df->Get(index,dfval);
                    if(val==0.0) {
                        eta+=  C  / h * dfval;
                    }
                    else {
                        eta-=  C  / h *dfval;
                    }
                    i++;
                }
                lastval = i -1;
            }
            else {
                eta=(*alpha)[i-1].value;
		    
                if(DEBUGLEVEL)
                    printf("Hit hinges: first %d last %d index %d\n",firstval,lastval,(*alpha)[i-1].index);
                for(int k=firstval;k<=lastval;k++) {
                    hinges.push_back((*alpha)[k].index);
                }
                
                break;
            }
        }
        if(i<=lastposindex) {
            // end because optimal eta found
            break;
        }
        else {
            // optimal eta not found, increase max step size
            max_stepsize = 2.0 * max_stepsize;
            if(DEBUGLEVEL)
                cout << "increase stepsize  to" <<  max_stepsize << "\n";
            if(lastval!=lastposindex) {
                // make sure that the  index / indices for the last eta
                // is stored whichas this  may be needed after the max step size is increased
                firstval =lastposindex ;
                lastval = lastposindex;
            }
            
            firstposindex = lastposindex;
            lastposindex =  numElements-1;
        }
    }
    //    Time4.Stop();
    
    return eta;

}



/**
   returns optimal step size starting from w in the
   direction of p. f=Xw with X beeing the data matrix.
   beta,df and alpha are temp variables of size of the
   data set which must be initialized by the caller
   Here normal sorting is used
*/
Scalar CSUBLBFGS::lineSearch(TheMatrix & p,TheMatrix & w,TheMatrix &f,TheMatrix * temp,TheMatrix * beta,TheMatrix * df,vector<ScalarWithIndex> * alpha,vector<int> &hinges)
{

    // 
    hinges.clear();
    //
    
    // start searching with interval that starts before
    // min_step and ends after min_step
    Scalar min_stepsize = -pow(10.0,-10.0); 
    
    lineSearchTime.Stop();
    lossAndGradientTime.Start();
    _data->ComputeF(p,*df);
    df->ElementWiseMult(_data->labels());
    lossAndGradientTime.Stop();
    lineSearchTime.Start();
    


    Scalar b;
    Scalar h;
    for (int i=0;i<numElements;i++) {
        Scalar val1;
        Scalar val2;
        f.Get(i,val1);
        df->Get(i,val2);
        (*alpha)[i].index = i;
        (*alpha)[i].value = (1- val1 ) / val2;
       
    }
    p.Dot(w,b);
    b*=lambda;
    
    p.Norm2(h);
    h=h*h;
    h*= lambda;
   
    // PRESORTING
    //put all elements smaller than min_stepsize before
    //the elements larger than min_stepsize
    //also stores largest number smaller than min_stepsize
    Scalar maxneg = 0.0; //   largest number smaller than min_stepsiz
    int maxnegindex = -1; // index of largest number smaller than min_stepsiz
    int firstposindex = -1; // all indices smaller than firstposindex are numbers smaller than min_stepsize
    
    // until we find a number larger than min_stepsize
    for(int i=0;i<numElements;i++) {
        Scalar thisval = (*alpha)[i].value;
        //int thisindex = (*alpha)[i].index;
        //store largest number
        if(thisval<min_stepsize) {
            if(thisval>maxneg || maxnegindex ==-1) {
                maxnegindex = i;
                maxneg = thisval;
            }
        }
        else {
            // first number found -> break
            firstposindex = i;
            break;

        }
    }
    // swapping will gives us the desired ordering
    for(int i=firstposindex+1;i<numElements;i++) {
        Scalar thisval = (*alpha)[i].value;
        
        if(thisval<min_stepsize) {
            if(thisval>maxneg || maxnegindex ==-1) {
                maxnegindex = i;
                maxneg = thisval;
            }
            // swap
            (*alpha)[i].swap((*alpha)[firstposindex]);
            if(maxneg==thisval) 
                maxnegindex = firstposindex;
            firstposindex++;
                
        }
    }


    if(maxnegindex!=-1) {
        // "move" largest number smaller then min_stepsize
        // one postion before all numbers larger then min_stepsize
         (*alpha)[maxnegindex].swap((*alpha)[firstposindex-1]);
   
    }
    
    // all numbers before firstposindex are irrelevant and
    // can have any kind of ordering

    if(firstposindex <0)
        firstposindex = 0;


    sort(alpha->begin() + firstposindex,alpha->end());
    
    

    int i= firstposindex + 1;
    Scalar eta = ((*alpha)[i-1].value +  (*alpha)[i].value) / 2.0;
    temp->Assign(f);
    temp->ScaleAdd(eta,*df);
    for (int j=0;j<numElements;j++) {
        Scalar val;
        temp->Get(j,val);
        val = (val<1);
        beta->Set(j,val);
    }
    if(DEBUGLEVEL)
        printf("starting index %d\n",i);
    //Time3.Stop();
   
    beta->Dot(*df,eta);
    eta = (C * eta - b) / h; 
    int firstval = 1;
    int lastval =1;
    
    while(i<numElements) {
     
        if(eta>=(*alpha)[i-1].value) {
            if(eta<=(*alpha)[i].value)
                break;
            Scalar flip = (*alpha)[i].value;
            firstval = i;
            while((*alpha)[i].value==flip) { // for all hinges with same alpha
                
                int index = (*alpha)[i].index;
                Scalar val;
                Scalar dfval;
                beta->Get(index,val);
                df->Get(index,dfval);
                if(val==0.0) {
                    eta+=  C / h * dfval;
                }
                else {
                    eta-=  C  / h *dfval;
                }
                i++;
            }
            lastval = i -1;
        }
        else {
            eta=(*alpha)[i-1].value;
            if(DEBUGLEVEL)
                printf("Hit hinges: first %d last %d index %d\n",firstval,lastval,(*alpha)[i-1].index);
            for(int k=firstval;k<=lastval;k++) {
                hinges.push_back((*alpha)[k].index);
            }
            
            break;
        }
    }

    
    return eta;

}





/* Find descent direction p given initial gradient g, parameters w,
   indices of the elements for which w is on the hinge (hinges)
   and the corresponding partial gradient (hingebetas).
   
   a2,gnew,temp1,temp2,ga are temp variables of size of the input
   space which must be initialized by the caller
*/
int CSUBLBFGS::getDescentDirection(TheMatrix & p,TheMatrix &a,TheMatrix & w,vector<int> &hinges,vector<double> &hingebetas,TheMatrix * a2, TheMatrix * gnew,TheMatrix * temp1,TheMatrix * temp2, TheMatrix * ga) {
   
    int checkduality = 1;


   
    const TheMatrix * labels = &_data->labels();
    
    int success = 1;

    gnew->Assign(a);

    
    Scalar gp;
    Scalar n1; 
    
   
    p.Assign(*gnew);
    BMult(p);
    p.Scale(-1.0);
    p.Norm2(n1);
    if (n1==0.0 || isnan(n1)) {
        printf("norm of p is invalid (%f)!\n",n1);
        success = 0;
    }
    gnew->Dot(p,gp);
    for(size_t i=0;i<hinges.size();i++) {
        // find gradient which leads to largest g'*p
        int index = hinges[i];
        Scalar val=0.0;
        Scalar label;
        // cerr << "index " << index;
        _data->ComputeFi(p,val,index);
        //p.Dot(*_data->x[index],val);
        labels->Get(index,label);
        val = C*val * label;
        
        // subtract "random" hinge values
        gp+=val*hingebetas[i];
        
        _data->AddElement(*gnew,index,C*label *hingebetas[i]);
        hingebetas[i] = 0.0;
        if(val<0) {
            gp-= val;
            _data->AddElement(*gnew,index,-label*C);
            hingebetas[i] = 1.0;
            }
     
    }
    if(DEBUGLEVEL)
        printf("gp %e \n",gp);
    Scalar gap=0.0;
    Scalar min_M=0.0;;
    //if(gp>0) {
       
        ga->Assign(a);
        if(DEBUGLEVEL)
            cout << "hinges size: " << hinges.size() << endl;

        if(checkduality) {
            // checking duality gap....
            
            ga->Dot(p,gap);
            min_M = gp - 0.5 * gap;
	    gap = min_M -  0.5 * gap;
        }
	// }
   
  
    
    int trynum = 1;
    
    while(gp>0 && trynum<MAXDESCENTTRY) { // not good direction
    //while(gp > 0 || gap > 1e-10 && trynum<MAXDESCENTTRY){
        if(DEBUGLEVEL)
            printf("trynum %d gp %e \n",trynum,gp);
        
        temp1->Assign(*ga);
        temp1->Minus(*gnew);
        temp2->Assign(*temp1);
        BMult(*temp2);
        Scalar nominator;
        Scalar denominator;
        Scalar mu;
        ga->Dot(*temp2,nominator);
        temp1->Dot(*temp2,denominator);
        mu = min((Scalar)1.0,nominator / denominator);

        if(mu<0) {
            printf("encountered negative mu %f  \n",mu);
            break;
        }
        ga->Scale(1.0-mu);
        ga->ScaleAdd(mu,*gnew);
        p.Scale(1.0-mu);
        temp2->Assign(*gnew);
        BMult(*temp2);
        temp2->Scale(-1.0);
        p.ScaleAdd(mu,*temp2);
        gnew->Dot(p,gp);
       

        for(size_t i=0;i<hinges.size() ;i++) {
            int index = hinges[i];
            Scalar val=0.0;
            Scalar label;
            _data->ComputeFi(p,val,index);
            
            labels->Get(index,label);
            val = val * label*C;
            if(val<0 && hingebetas[i]==0) {

                gp-= val;
                _data->AddElement(*gnew,index,-label*C);
                hingebetas[i]=1.0;
            }
            if(val>0 && hingebetas[i]==1) {
                 
                gp+= val;
                //gnew->ScaleAdd(label / ((double) numElements),*_data->x[index]);
                _data->AddElement(*gnew,index,label*C);
                hingebetas[i]=0;
             }
          
        }
        temp1->Assign(p);
        BMult(*temp1);
        Scalar pnorm;
        temp1->Dot(p,pnorm);


        // Duality Gap
        if(checkduality) {
            ga->Dot(p,gap);
            if( gp - 0.5*gap < min_M)
                min_M =  gp - 0.5*gap;
            gap = min_M - 0.5*gap;
	}
        if(DEBUGLEVEL)
            cout << "try " << trynum  << "gp " << gp <<" pnorm " << pnorm  << "ob " << (pnorm * 0.5 + gp) << " mu "  << mu << " gap "  << gap << "\n" ;
        trynum++;
    } 
    if(trynum==MAXDESCENTTRY) 
        success= 0;
    if(USENEWGRADIENT && trynum > 1)
        a.Assign(*ga);


    if(trynum>1)
        printf("more than 1 inner loop executed : try %d    \n",trynum);
   
    return success;
}


void CSUBLBFGS::Train()
{
   
    TheMatrix &w = _model->GetW();
    
    int dimOfW = _data->dim();
  
    b_start =0;
    actual_bsize=0;
    int    iter = 0;              // iteration count
    Scalar old_obj;
    Scalar obj = 0.0;            // objective function value        
    Scalar loss = 0.0;            // loss function value        
    //int pos =0;
   
    vector<int> hinges;
    vector<double> hingebetas;
    
    
    // line search tmps
    TheMatrix * temp = new TheMatrix(numElements,1);
    TheMatrix * beta2 = new TheMatrix(numElements,1);
    TheMatrix * df = new TheMatrix(numElements,1);
    vector<ScalarWithIndex> alpha(numElements);
    
    // descent direction tmps
    TheMatrix * a3 = new TheMatrix(dimOfW,1);
    
    TheMatrix * gnew= new TheMatrix(dimOfW,1);
    TheMatrix * temp1= new TheMatrix(dimOfW,1);
    TheMatrix * temp2= new TheMatrix(dimOfW,1);
    TheMatrix * ga= new TheMatrix(dimOfW,1);
    


   TheMatrix* a = 0;                // gradient
   TheMatrix * p =0;
   TheMatrix* a2 = 0;
   TheMatrix * f = 0;
   TheMatrix * beta = 0;
   
  
   // check and confirm program parameters from CInfo and/or Configuration file
   ConfirmProgramParameters();



   AllocateBuffer(dimOfW);
   a = new TheMatrix(dimOfW, 1);
   a2 = new TheMatrix(dimOfW, 1);
   p = new TheMatrix(dimOfW,1);
   f = new TheMatrix(numElements,1);
   beta = new TheMatrix(numElements,1);
   

   // start training
   totalTime.Start();
   
 
   // gradient calculation
   lossAndGradientTime.Start();
   evaluateFunction(w, loss,obj,*a,hinges,hingebetas,*beta,*f,1);
   lossAndGradientTime.Stop();
   double eta = 0.0;
   while(1) {
       
       iter++;
      
    
       descentDirectionTime.Start();
       int success = getDescentDirection(*p,*a,w, hinges,hingebetas,a3, gnew,temp1, temp2,ga);
       descentDirectionTime.Stop();
       
     
       if(!success) {
           printf("no descent direction found, aborting\n"); 
           printf("it %d, obj %.20e, loss %e, eta %e \n",iter,obj,loss,eta);
           break;
       }
           
   
       lineSearchTime.Start();
       eta =lineSearch2(*p,w,*f,temp,beta2,df,& alpha,hinges);
       lineSearchTime.Stop();
       if(eta==0) {
           printf("encountered eta 0, aborting \n");
	   printf("it %d, obj %.20e, loss %e, eta %e \n",iter,obj,loss,eta);
           break;
       }
      old_obj = obj; 
           
      w.ScaleAdd(eta,*p);
      f->ScaleAdd(eta,*df);
  

      if(0) {
          // something needed against numerical problems
          for(size_t i=0;i<hinges.size();i++) {
              int index =hinges[i];
              f->Set(index,1.0);
          }
      }
      lossAndGradientTime.Start();
      
      evaluateFunction(w, loss,obj,*a2,hinges,hingebetas,*beta,*f,0);
 
      lossAndGradientTime.Stop();
   
      
      p->Scale(eta);
      a->Scale(-1.0);
      a->Add(*a2);

      totalTime.Stop();
      //printf("eta %e it %d, obj %.20e, loss %e, time %e\n",eta,iter,obj,loss,totalTime.CPUTotal()); 
      printf("it %d, obj %.20e, time %e, loss %e, eta %e\n",iter,obj,totalTime.CPUTotal(), loss, eta); 
      totalTime.Start();
      fflush(stdout);
      /* buffer management */
      AddToBuffer(*p,*a);

      a->Assign(*a2);
      
      if(iter >= maxNumOfIter)
      { 
         printf("\nWARNING: program exceeded maximum number of iterations (%d) !\n", maxNumOfIter);
         //printf("eta %e it %d, obj %.20e, loss %e\n",eta,iter,obj,loss);
	 printf("it %d, obj %.20e, loss %e, eta %e \n",iter,obj,loss,eta);
         break;
      }
   }

   printf("\nLegends: \neta: optimal step size \nit: iteration number \nobj: objective value\nloss loss function value\ntime: CPU time used\n\n"); 
   
   totalTime.Stop();
   // end of training

   // display timing profile
   printf("\nCPU seconds in:\n");
   printf("1. loss and gradient : %8.2f\n", lossAndGradientTime.CPUTotal());
   printf("2. line search       : %8.2f\n", lineSearchTime.CPUTotal());
   printf("3. descent direction : %8.2f\n", descentDirectionTime.CPUTotal());
   printf("                Total: %8.2f\n", totalTime.CPUTotal());

  
   delete a;              
   delete p;
   delete a2;
   delete f;
   delete beta;
   
   delete temp;
   delete beta2;
   delete df;
   

   delete a3;
   delete gnew;
   delete temp1;
   delete temp2;
   delete ga;
   
}


#endif
