/* Copyright (c) 2006, National ICT Australia 
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
 * Authors: S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
 *
 * Created: (17/01/2008) 
 *
 * Last Updated:
 */

#ifndef _FACELOSS_CPP_
#define _FACELOSS_CPP_

#include "faceloss.hpp"
#include "sml.hpp"

CFaceLoss::CFaceLoss(CModel* &model, CVecData* &data)
  :CScalarLoss(model, data) {
  
  // do nothing serious here at the moment...                       
  std::cout << "In FaceLoss!" << std::endl;
  return;
}

        
/**  
 * Display (usage) message.
 */
void CFaceLoss::Usage()
{
  std::cout << "FaceLoss: no parameters" << std::endl; 
}


/**  
 * Do nothing for now.
 */
void CFaceLoss::Loss(Scalar& loss, TheMatrix& f)
{
  return;
}

/**  
 *  Compute loss and partial derivative of Face loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CFaceLoss::LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l){
  
  // svnvish: BUGBUG
  // Ignore Y for now.
  // Stupid assumption: The first face of an object is the true face. 
  
  loss = 0.0;	
  l.Zero();  
  int num_obj = _data->NumOfSubset();
  Scalar* f_array = f.Data();  
  Scalar* y_array = _data->labels().Data();
  Scalar* l_array = l.Data();
  for(int i=0; i<num_obj; i++){
    int offset = _data->subset[i].startIndex;
    int subsetsize = _data->subset[i].size;
    // At least two faces per object.
    assert(subsetsize>1);
    Scalar max_fx = f_array[1+offset];
    int max_idx = 1+offset;
    for(int i = 1; i < subsetsize; i++){
      if(f_array[i+offset] > max_fx){
        max_fx = f_array[i+offset];
        max_idx = i+offset;
      }
    }

    // enforce a margin of 1
    if(f_array[offset] < 1+max_fx){
      // We incur a loss so we need to set loss and gradient
      loss +=  (1+ max_fx - f_array[offset]);
      l_array[offset] = -1;
      l_array[max_idx] = 1;
    }
    
  }
}

/** 
 * Transform function value f := <w,x> into label
 *
 *  @param f [read/write] function values / labels
 */
void CFaceLoss::Transform(TheMatrix& f)
{
  
	Scalar* f_array = f.Data(); 
  int num_obj = _data->NumOfSubset();
  for(int i=0; i<num_obj; i++){
    int offset = _data->subset[i].startIndex;
    int subsetsize = _data->subset[i].size;
    Scalar max_fx = f_array[offset];
    int max_idx = 0;
    for(int i = 0; i < subsetsize; i++){
      if(f_array[i+offset] > max_fx){
        max_fx = f_array[i+offset];
        max_idx = i;
      }
    }
    for(int i = 0; i < subsetsize; i++){
      f_array[i+offset] = max_idx;
    }
  }
  return;
}

/** Evaluate the performance of the model on _data
 *
 *  @param f [read] function values
 *  @param predict [write] prdicted labels
 */
void CFaceLoss::Perf(TheMatrix& f, TheMatrix& predict)
{
  Scalar* f_array = f.Data(); 
  Scalar* p_array = predict.Data(); 
  int num_obj = _data->NumOfSubset();
  for(int i=0; i<num_obj; i++){
    int offset = _data->subset[i].startIndex;
    int subsetsize = _data->subset[i].size;
    Scalar max_fx = f_array[offset];
    int max_idx = 0;
    for(int i = 0; i < subsetsize; i++){
      if(f_array[i+offset] > max_fx){
        max_fx = f_array[i+offset];
        max_idx = i;
      }
    }
    for(int i = 0; i < subsetsize; i++){
      p_array[i+offset] = max_idx;
    }
  }
  return;

  return;
}

#endif
