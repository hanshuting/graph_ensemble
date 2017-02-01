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
 * Created: (02/11/2007) 
 *
 * Last Updated: (13/11/2007)   
 */

#ifndef _SOFTMARGINLOSS_HPP_
#define _SOFTMARGINLOSS_HPP_

#include "sml.hpp"
#include "binaryclassificationloss.hpp"
#include "model.hpp"

/**  
 * Class to represent binary soft margin loss function:
 * loss = max(0, margin - y*f(x))
 * where f(x) := <w, x> and margin \geq 0. By default margin = 1.
 */
class CSoftMarginLoss : public CBinaryClassificationLoss 
{
protected:
	void Loss(Scalar& loss, TheMatrix& f);
	void LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l);
	
public:    
	
	CSoftMarginLoss(CModel* &model, CVecData* &data)
		: CBinaryClassificationLoss(model, data) {}
	virtual ~CSoftMarginLoss() {}
	virtual void Usage();
      
};

#endif
