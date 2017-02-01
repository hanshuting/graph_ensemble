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

#ifndef _NDCGRANKLOSS_HPP_
#define _NDCGRANKLOSS_HPP_

#include "sml.hpp"
#include "rankloss.hpp"
#include "model.hpp"

/** NDCG Ranking loss.
 *  Will create a superclass such as CRankLoss for all ranking loss soon...
 */
class CNDCGRankLoss : public CRankLoss
{
protected:               
  int truncation;
  vector<int> current_ideal_pi;
  vector< vector<int> > sort_vectors; /* sort indices per query */
  int max_subset_size;
  vector<double> c;
  int *pi;
  vector < vector<double> > bs;
  vector<double> a;

  void Loss(Scalar& loss, TheMatrix& f);
  void LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l);
  void delta(int size, vector<double> a, vector<double> b, 
	     int *pi, double &value);
  void find_permutation(int size, int offset, vector<double> a, 
			vector<double> b, vector<double> c, 
			Scalar *f, int *pi);
  void compute_coefficients(int offset, int size, Scalar *y, 
			    vector<int> ideal_pi, 
			    vector<double> a, vector<double> b);
  double get(Scalar *f, int offset, int i);
  void add(TheMatrix &l, int offset, int i, double value);
  void compute_b(int offset, int size, Scalar *y, 
		 vector<int> ideal_pi, 
		 vector<double> &input_b);
  void compute_a(int size, vector<double> &input_a);

public:    
	
	CNDCGRankLoss(CModel* &model, CVecData* &data);
	
	virtual ~CNDCGRankLoss() {if (pi){free(pi);pi = 0;}}
	virtual void Usage();
      
};

#endif
