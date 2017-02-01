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

#ifndef _ORDINALREGRESSIONLOSS_HPP_
#define _ORDINALREGRESSIONLOSS_HPP_

#include "sml.hpp"
#include "rankloss.hpp"
#include "model.hpp"

/** OrdinalRegression loss.
 *  Will create a superclass such as CRankLoss for all ranking loss soon...
 */
class COrdinalRegressionLoss : public CRankLoss 
{
protected:               
  vector<int> current;
  vector< vector<int> > sort_vectors; /* sort indices per query */
	void Loss(Scalar& loss, TheMatrix& f);
	void LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l);
  
public:    
	
	COrdinalRegressionLoss(CModel* &model, CVecData* &data);
	
	virtual ~COrdinalRegressionLoss() {}
	virtual void Usage();
      
};

#endif
