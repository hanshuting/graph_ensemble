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
 * Last Updated: (19/11/2007)   
 */

#include "epsiloninsensitiveloss.hpp"
#include "configuration.hpp"
#include "sml.hpp"

/** Constructor
 *
 *  @param w [read] weight vector
 *  @param data [read] pointer to dataset
 */
CEpsilonInsensitiveLoss::CEpsilonInsensitiveLoss(CModel* &model, CVecData* &data)
	: CUnivariateRegressionLoss(model, data),
	  epsilon(0.1)
{
	// read loss function parameters
	Configuration &config = Configuration::GetInstance();
	
	if(config.IsSet("EpsilonInsensitive.epsilon")) {
		epsilon = config.GetDouble("EpsilonInsensitive.epsilon");
	}
}

/**  
 * Display (usage) message
 */
void CEpsilonInsensitiveLoss::Usage()
{
  std::cout << "Epsilon Insensitive Regression Loss parameters:" << std::endl; 
  std::cout << "   Type:   Key:                        Range:      Default:" << std::endl;
  std::cout << "   double  EpsilonInsensitive.epsilon  [0,\\infty)  0.1" << std::endl; 
}



/**  
 *  Compute epsilon insensitive regression loss. CAUTION: f is passed by
 *  reference and is changed within this function. This is done for
 *  efficiency reasons, otherwise we would have had to create a new copy
 *  of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CEpsilonInsensitiveLoss::Loss(Scalar& loss, TheMatrix& f)
{
	loss = 0;
	f.Minus(_data->labels()); // f = f - y
	Scalar* f_array = f.Data(); 
	int len = f.Length();
	
	for(int i=0; i < len; i++)
		loss += std::max((Scalar)0.0, SML::abs(f_array[i]) - epsilon);
}

/**  
 *  Compute loss and gradient of epsilon insensitive regression
 *  loss. CAUTION: f is passed by reference and is changed within this
 *  function. This is done for efficiency reasons, otherwise we would
 *  have had to create a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 *  @param l [write] partial derivative of loss function w.r.t. f
 */
void CEpsilonInsensitiveLoss::LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l)
{
	f.Minus(_data->labels()); // f = f - y
	Scalar* f_array = f.Data();  // pointer to memory location of f (faster element access)
	int len = f.Length();
	
	l.Zero();  // grad := l'*X
	
	Scalar loss_i = 0;
	for(int i=0; i < len; i++) 
	{
		loss_i = std::max((Scalar)0.0, SML::abs(f_array[i]) - epsilon);
		loss += loss_i;
		if(loss_i > 0.0)
			l.Set(i, SML::sgn(f_array[i]));
	}
}
