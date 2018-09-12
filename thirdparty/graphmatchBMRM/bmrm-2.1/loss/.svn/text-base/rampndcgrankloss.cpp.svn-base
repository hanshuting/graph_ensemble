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

#ifndef _RAMPNDCGRANKLOSS_CPP_
#define _RAMPNDCGRANKLOSS_CPP_

#include "rampndcgrankloss.hpp"
#include "sml.hpp"
#include <ext/numeric>         // for iota

CRampNDCGRankLoss::CRampNDCGRankLoss(CModel* &model, CVecData* &data)
  : CRankLoss(model, data) 
{
  truncation = 10;
  // do nothing serious here at the moment...                       
  std::cout << "In CRampNDCGRankLoss!" << std::endl;
  int num_query = _data->NumOfSubset();
  Scalar *Y_array = _data->labels().Data();
  max_subset_size = 0;
  for (int q=0;q<num_query;q++){
    int offset = _data->subset[q].startIndex;
    int subsetsize = _data->subset[q].size;
    if (max_subset_size < subsetsize){
      max_subset_size = subsetsize;
    }
    int length = subsetsize;
    vector<int> indices(length);
    iota(indices.begin(), indices.end(), 0);
    indirect_greater_than<Scalar> igt(&(Y_array[offset]));
    //indirect_less_than<Scalar> igt(&(Y_array[offset]));
    sort(indices.begin(), indices.end(),igt);
    sort_vectors.push_back(indices);

  }
  // use largest query information to compute c_i
  for (int i=0;i<max_subset_size;i++){
    c.push_back( max(truncation + 1 - i,0) );
  }


}

        
/**  
 * Display (usage) message.
 */
void CRampNDCGRankLoss::Usage()
{
  std::cout << "RampNDCGRankLoss: parameters description will be in soon..." << std::endl; 
}


void CRampNDCGRankLoss::find_permutation(int size, int offset, 
				     Scalar *f, vector<int> &input_pi){

  // for this ramp loss, permutation is the sorting order
  iota(input_pi.begin(), input_pi.end(), 0);
  indirect_greater_than<Scalar> igt(&(f[offset]));
  //indirect_less_than<Scalar> igt(&(f[offset]));
  sort(input_pi.begin(), input_pi.end(),igt);

}



/**  
 *  Compute RampNDCGRank loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CRampNDCGRankLoss::Loss(Scalar& loss, TheMatrix& f)
{
  // chteo: here we make use of the subset information 
        
  loss = 0.0;	
  Scalar* f_array = f.Data();  
  for(int q=0; q < _data->NumOfSubset(); q++)
    {
      int offset = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      current_ideal_pi = sort_vectors[q];
      /* find the best permutation */
      vector<int> pi(subsetsize);
      find_permutation(subsetsize, offset, f_array, pi);
      
      /* compute the loss */
      for (int i=0;i<subsetsize;i++){
	loss = loss + c[i]*(get(f_array, offset, pi[i]) - get(f_array, offset, i));
      }
      
    }

}

/**  
 *  Compute loss and partial derivative of RampNDCGRank loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CRampNDCGRankLoss::LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l)
{
  // chteo: here we make use of the subset information 
        
  loss = 0.0;	
  l.Zero();  
  Scalar* f_array = f.Data();  
  for(int q=0; q < _data->NumOfSubset(); q++)
    {
      int offset = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      current_ideal_pi = sort_vectors[q];
      /* find the best permutation */
      vector<int> pi(subsetsize);
      find_permutation(subsetsize, offset, f_array, pi);
      for (int i=0;i<subsetsize;i++){
	loss = loss + c[i]*(get(f_array, offset, pi[i]) - get(f_array, offset, i));
      }
      
      for (int i=0;i<subsetsize;i++){
	add(l, offset, i, -c[i]);
	add(l, offset, pi[i], +c[i]);
      }
    }
  

}

double CRampNDCGRankLoss::get(Scalar *f, int offset, int i){
  return f[offset + current_ideal_pi[i]];
}

void CRampNDCGRankLoss::add(TheMatrix &l, int offset, int i, double value){
  Scalar temp;
  l.Get(offset + current_ideal_pi[i], temp);
  l.Set(offset + current_ideal_pi[i], temp + value);
}

#endif
