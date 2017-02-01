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
 * Authors: Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (06/01/2008) 
 *
 * Last Updated:
 */

#ifndef _PREFERENCERANKLOSS_CPP_
#define _PREFERENCERANKLOSS_CPP_

#include "preferencerankloss.hpp"
#include "sml.hpp"

CPreferenceRankLoss::CPreferenceRankLoss(CModel* &model, CVecData* &data)
        : CRankLoss(model, data) 
{
       // do nothing serious here at the moment...                       
       std::cout << "In CPreferenceRankLoss!" << std::endl;
}

        
/**  
 * Display (usage) message.
 */
void CPreferenceRankLoss::Usage()
{
        std::cout << "PreferenceRankLoss: parameters description will be in soon..." << std::endl; 
}



/**  
 *  Compute PreferenceRank loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CPreferenceRankLoss::Loss(Scalar& loss, TheMatrix& f)
{
  const TheMatrix &Y = _data->labels();
  Scalar* f_array = f.Data(); 
  for(int q=0; q < _data->NumOfSubset(); q++)
    {
      int k = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      for(int i=k; i < k+subsetsize; i++)
	{
	  for(int j=k; j < k+subsetsize; j++){
	    Scalar yi;
	    Scalar yj; 
	    Y.Get(i, yi);
	    Y.Get(j, yj);
	    if (yi < yj){
	      if (f_array[j] - f_array[i] < 1){
		loss += (yj - yi)*(1 + f_array[i] - f_array[j]);
	      }
	    }
	  }
	}
    }
 
}

/**  
 *  Compute loss and partial derivative of PreferenceRank loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CPreferenceRankLoss::LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l)
{
  const TheMatrix &Y = _data->labels();
  Scalar* f_array = f.Data(); 
  l.Zero();  
  
  for(int q=0; q < _data->NumOfSubset(); q++)
    {
      int k = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      for(int i=k; i < k+subsetsize; i++)
	{
	  for(int j=k; j < k+subsetsize; j++){
	    Scalar yi;
	    Scalar yj; 
	    Y.Get(i, yi);
	    Y.Get(j, yj);
	    if (yi < yj){
	      if (f_array[j] - f_array[i] < 1){
		loss += (yj - yi)*(1 + f_array[i] - f_array[j]);
		Scalar temp; 
		l.Get(i,temp);
		l.Set(i,temp + (yj - yi));
		l.Get(j,temp);
		l.Set(j, temp - (yj - yi));
	      }
	    }
	  }
	}
    }
}

#endif
