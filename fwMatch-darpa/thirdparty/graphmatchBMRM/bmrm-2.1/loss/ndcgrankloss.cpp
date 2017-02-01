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

#ifndef _NDCGRANKLOSS_CPP_
#define _NDCGRANKLOSS_CPP_

#include "ndcgrankloss.hpp"
#include "sml.hpp"
#include "lap.hpp"
#include "numeric_wrap.h"         // for iota

CNDCGRankLoss::CNDCGRankLoss(CModel* &model, CVecData* &data)
  : CRankLoss(model, data)
{
  truncation = 10;
  // do nothing serious here at the moment...
  std::cout << "In CNDCGRankLoss!" << std::endl;
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

    // use this to compute bs
    vector<double> b;
    compute_b(offset, subsetsize, Y_array, indices, b);
    bs.push_back( b );
  }
  // use largest query information to compute c_i
  for (int i=0;i<max_subset_size;i++){
    //c.push_back( 1.0/pow((i+1.0),1.4) );
    c.push_back( max(truncation+1-i,0) );
  }
  pi = (int *) malloc(sizeof(int)*max_subset_size);
  compute_a(max_subset_size, a);
  //cout << "finished init\n";
  //cout << "max_subset_size " << max_subset_size <<" \n";
  //cout << "a.size " << a.size() <<" \n";
  //cout << "c.size " << c.size() <<" \n";

}


/**
 * Display (usage) message.
 */
void CNDCGRankLoss::Usage()
{
  std::cout << "NDCGRankLoss: parameters description will be in soon..." << std::endl;
}



void CNDCGRankLoss::compute_coefficients(int offset, int size, Scalar *y,
					 vector<int> ideal_pi,
					 vector<double> a, vector<double> b){

  double dy;
  eval_dcg(offset, size, truncation, y, ideal_pi, dy);
  for (int i=0;i<size;i++){
    if (i<truncation){
      // svnvish: BUGBUG
      // does not compile on my machine
      // a[i] = 1.0/log2(i+2.0);
    }
    else{
      a[i] = 0.0;
    }
    if (dy == 0){
      b[i] = 1;
    }else{
      b[i] = (pow(2.0,(double) y[offset + ideal_pi[i]]) - 1.0)/dy;
    }
  }
}

void CNDCGRankLoss::compute_b(int offset, int size, Scalar *y,
			      vector<int> ideal_pi,
			      vector<double> &input_b){

  double dy;
  eval_dcg(offset, size, truncation, y, ideal_pi, dy);
  for (int i=0;i<size;i++){
    if (dy == 0){
      input_b.push_back(1);
    }else{
      input_b.push_back( (pow(2.0,(double) y[offset + ideal_pi[i]]) - 1.0)/dy );
    }
  }
}

void CNDCGRankLoss::compute_a(int size, vector<double> &input_a){
  //cout <<" size = "<< size << "\n";
  for (int i=0;i<size;i++){
    if (i<truncation){
      // svnvish: BUGBUG
      // does not compile on my machine
      // input_a.push_back( 1.0/log2(i+2.0) );
    }
    else{
      input_a.push_back(0.0);
    }
  }
  //cout << "input-a.size after compute_a "<< input_a.size() << "\n";
}

void CNDCGRankLoss::find_permutation(int size, int offset, vector<double> a,
				     vector<double> b, vector<double> c,
				     Scalar *f, int *input_pi){
  /* setting up C_ij */
  double **C = (double **)malloc(sizeof(double*)*size);
  for (int i=0;i<size;i++){
    C[i] = (double *)malloc(sizeof(double)*size);
  }
  //cout << "before inner loop\n";
  //cout << "size = "<< size << "\n";
  //cout << "a.size = "<< a.size() <<"\n";
  //cout << "b.size = "<< b.size() << "\n";
  //cout << "c.size = "<< c.size() << "\n";
  for (int i=0;i<size;i++){
    for (int j=0;j<size;j++){
      C[i][j] = -c[i]*get(f, offset, j) + b[j]*a[i];
    }
  }
  //cout << "after inner loop\n";
  int *row = (int *)malloc(sizeof(int)*size);

  //cout << "before lap\n";
  lap(size,C,input_pi,row);
  //cout << "after lap\n";

  free(row);
  for (int i=0;i<size;i++){
    free(C[i]);
  }
  free(C);
}

void CNDCGRankLoss::delta(int size, vector<double> a, vector<double> b,
			  int *pi, double &value){
  value = 0;
  for (int i=0;i<size;i++){
    value += (b[i]- b[pi[i]])*a[i];
    //value += - b[i]*a[pi[i]];
  }
  if (value < 0){
    //cout << "\nvalue = "<< value<< "\n";
  }
  else if (value > 1){
    //cout << "\nvalue = "<< value<< "\n";
  }
}


/**
 *  Compute NDCGRank loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f.
 *
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector.
 */
void CNDCGRankLoss::Loss(Scalar& loss, TheMatrix& f)
{
  // chteo: here we make use of the subset information

  loss = 0.0;
  Scalar* f_array = f.Data();
  for(int q=0; q < _data->NumOfSubset(); q++)
    {
      int offset = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      current_ideal_pi = sort_vectors[q];
      vector<double> b = bs[q];

      //compute_coefficients(offset, subsetsize, y_array, current_ideal_pi, a, b);

      /* find the best permutation */
      find_permutation(subsetsize, offset, a, b, c, f_array, pi);

      /* compute the loss */
      double value;
      delta(subsetsize, a, b, pi, value);

      loss += value;

      for (int i=0;i<subsetsize;i++){
	loss = loss + c[i]*(get(f_array, offset, pi[i]) - get(f_array, offset, i));
      }
      //free(c);
      //free(a);
      //free(b);
      //free(pi);

    }

}

/**
 *  Compute loss and partial derivative of NDCGRank loss w.r.t f
 *
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CNDCGRankLoss::LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l)
{
  // chteo: here we make use of the subset information

  loss = 0.0;
  l.Zero();
  Scalar* f_array = f.Data();
  for(int q=0; q < _data->NumOfSubset(); q++)
    {
      //cout << "q = "<< q <<endl;
      int offset = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      current_ideal_pi = sort_vectors[q];
      vector<double> b = bs[q];

      //compute_coefficients(offset, subsetsize, y_array, current_ideal_pi, a, b);

      //cout << "before finding permutation\n";
      /* find the best permutation */
      find_permutation(subsetsize, offset, a, b, c, f_array, pi);
      //cout << "after finding permutation\n";

      //cout << "before finding delta\n";
      /* compute the loss */
      double value;
      delta(subsetsize, a, b, pi, value);
      //cout << "before finding delta\n";

      loss += value;

      for (int i=0;i<subsetsize;i++){
	loss = loss + c[i]*(get(f_array, offset, pi[i]) - get(f_array, offset, i));
      }

      for (int i=0;i<subsetsize;i++){
	//add(l, offset, i, c[pi[i]] - c[i]);
	add(l, offset, i, - c[i]);
	add(l, offset, pi[i], c[i]);
      }
    }


}

double CNDCGRankLoss::get(Scalar *f, int offset, int i){
  return f[offset + current_ideal_pi[i]];
}

void CNDCGRankLoss::add(TheMatrix &l, int offset, int i, double value){
  Scalar temp;
  l.Get(offset + current_ideal_pi[i], temp);
  l.Set(offset + current_ideal_pi[i], temp + value);
}

#endif
