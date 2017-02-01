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

#ifndef _FACERANKLOSS_HPP_
#define _FACERANKLOSS_HPP_

#include "sml.hpp"
#include "scalarloss.hpp"
#include "model.hpp"

/** 
 * Quick and dirty hack to test if we can predict the correct face on
 * which an object rests on.
 */
class CFaceLoss : public CScalarLoss
{
protected:               
	void Loss(Scalar& loss, TheMatrix& f);
	void LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l);
  void Transform(TheMatrix& f);
  void Perf(TheMatrix& f, TheMatrix& predict);
	
public:    
	
	CFaceLoss(CModel* &model, CVecData* &data);
	
	virtual ~CFaceLoss() {}
	virtual void Usage();
  
};

#endif
