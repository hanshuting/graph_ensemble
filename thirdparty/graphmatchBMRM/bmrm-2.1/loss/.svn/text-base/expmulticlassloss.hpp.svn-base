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
 *          Choon Hui Teo (Choonhui.Teo@anu.edu.au)
 *
 * Created: (28/11/2007) 
 *
 * Last Updated: (06/011/2008)   
 */

#ifndef _EXPMULTICLASSLOSS_HPP_
#define EXPMULTICLASSLOSS_HPP_

#include <string>
#include "sml.hpp"
#include "data.hpp"
#include "vecdata.hpp"
#include "loss.hpp"
#include "model.hpp"
#include "timer.hpp"

/**  
 * Class for encapsulating the following ranking loss:
 *
 * loss = g(w|x) - <\phi(x,y), w>
 * 
 * where \phi(x, y) := [ \phi_a(x).y, \phi_b(x).y^2]^T
 * 
 * and g(w|x) = log \sum_{y} exp(<\phi(x,y), w>)
 * 
 * Of course the assumption is that w := [w_a, w_b]
 * 
 * NOTE: We assume that the labels are integers which range from 1 to
 * the data.maxLabel (typically assumed to be small, of the order of 5
 * or so). 
 * 
 */
class CExpMulticlassLoss : public CLoss 
{               
protected:
        Scalar* score;
        Scalar* mygrad;
        Scalar* val;
        TheMatrix* f;
        TheMatrix* l;
        TheMatrix* mw;
        TheMatrix* mgrad;
        
	/** pointer to dataset
	 */
	CVecData* _data;
    
	/** shorthand for number of examples
	 */
	unsigned int m;
    
	/** number of classes
	 */
	unsigned int classes;
    
	/** number of hyperplane
	 */
	unsigned int numOfW;
    
	/** dimensionality of each hyperplane
	 */
	unsigned int dimOfW;
	
	/**  Timer to keep track of the time spent on computing gradient
         */
        CTimer gradtime;

	
	// Internal helper functions
    
	/** Calculate Wx
	 */
	void CalcF(const TheMatrix& w, TheMatrix& f);
    
	/** Turn Wx into a probability distribution
	 */
	void CalcProb(const Scalar* fval, Scalar* prob);
	
public:
        CExpMulticlassLoss() : score(0), mygrad(0), val(0), f(0), l(0), mw(0), mgrad(0), _data(0), m(0), classes(0), numOfW(2), dimOfW(0) {} 
	CExpMulticlassLoss(CModel* &model, CVecData* &data); 
	virtual ~CExpMulticlassLoss();
        
	// Methods
	virtual void Usage();
	virtual void ComputeLoss(Scalar& loss);
	virtual void ComputeLossAndGradient(Scalar& loss, TheMatrix& grad);
	virtual void Predict(CModel* model);
	virtual void Evaluate(CModel* model);
    
};

#endif
