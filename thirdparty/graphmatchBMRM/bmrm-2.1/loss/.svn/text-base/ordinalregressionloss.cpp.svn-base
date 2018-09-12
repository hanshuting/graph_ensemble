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

#ifndef _ORDINALREGRESSIONLOSS_CPP_
#define _ORDINALREGRESSIONLOSS_CPP_

#include "ordinalregressionloss.hpp"
#include "sml.hpp"
#include <ext/numeric>         // for iota

COrdinalRegressionLoss::COrdinalRegressionLoss(CModel* &model, CVecData* &data)
        : CRankLoss(model, data) 
{
       // do nothing serious here at the moment...                       
       std::cout << "In COrdinalRegressionLoss!" << std::endl;

       int num_query = _data->NumOfSubset();
       Scalar *Y_array = _data->labels().Data();
       for (int q=0;q<num_query;q++){
	 int k = _data->subset[q].startIndex;
	 int subsetsize = _data->subset[q].size;
	 int length = subsetsize;
	 vector<int> indices(length);
	 iota(indices.begin(), indices.end(), 0);
	 indirect_greater_than<Scalar> igt(&(Y_array[k]));
	 sort(indices.begin(), indices.end(),igt);
	 sort_vectors.push_back(indices);
       }
       
}

        
/**  
 * Display (usage) message.
 */
void COrdinalRegressionLoss::Usage()
{
        std::cout << "OrdinalRegressionLoss: parameters description will be in soon..." << std::endl; 
}



/**  
 *  Compute OrdinalRegression loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void COrdinalRegressionLoss::Loss(Scalar& loss, TheMatrix& f)
{
  //const TheMatrix &Y = _data->labels();
  //Scalar* f_array = f.Data(); 
 
}

/**  
 *  Compute loss and partial derivative of OrdinalRegression loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void COrdinalRegressionLoss::LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l)
{
  const TheMatrix &Y = _data->labels();
  Scalar* f_array = f.Data(); 
  l.Zero();  
  loss = 0;
  for(int q=0; q < _data->NumOfSubset(); q++)
    {
      int offset = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      current = sort_vectors[q];

      int maxy = 0;
      for (int i=0;i<subsetsize;i++){
	Scalar yi;

	Y.Get(offset + i, yi);
	if (yi>maxy){
	  maxy = (int)yi;
	}
      }
      int labelsize = maxy + 1;
      vector<int> U(labelsize);
      vector<int> L(labelsize);
      for (int i=0;i<labelsize;i++){
	L[i] = 0;
	U[i] = 0;
      }
      for (int j=0;j<subsetsize;j++){
	for (int i=0;i<labelsize;i++){
	  Scalar yj;
	  Y.Get(offset + j, yj);
	  int myy = (int)yj;
	  if (myy==i){
	    U[i] += 1;
	    break;
	  }
	}
      }
      double M = subsetsize*subsetsize;
      for (int i=0;i<labelsize;i++){
	M -= U[i]*U[i];
      }
      M = 0.5*M;
      if (M != 0){
	/* initialization of loss and gradient vector */
	vector<int> indices(2*subsetsize);
	vector<Scalar> c(2*subsetsize);
	for (int i=0;i<subsetsize;i++){
	  c[i] = f_array[offset + current[i]] - 0.5;
	  c[i + subsetsize] = f_array[offset + current[i]] + 0.5;
	}
	/* sort the indices with descending order */
	indirect_less_than<Scalar> ilt(&(c[0]));
	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(), ilt);
	
	for (int i=0;i<2*subsetsize;i++){
	  int j = indices[i]%subsetsize;
	  Scalar yj;
	  Y.Get(offset + current[j], yj);
	  int z = (int)yj;
	  
	  if (indices[i] < subsetsize){
	    for (int k=0;k<z;k++){
	      loss = loss - (z - k)*U[k]*c[j]/M;
	      Scalar temp;
	      l.Get(offset + current[j], temp);
	      l.Set(offset + current[j], temp - (z - k)*U[k]/M);
	    }
	    L[z] = L[z] + 1;
	  }
	  else{
	    for (int k=z+1;k<labelsize;k++){
	      loss = loss + (k - z)*L[k]*c[j+m]/M;
	      Scalar temp;
	      l.Get(offset + current[j], temp);
	      l.Set(offset + current[j], temp + (k - z)*L[k]/M);
	    }
	    U[z] = U[z] - 1; 
	  }
	}
	
      }

    }
}


#endif
