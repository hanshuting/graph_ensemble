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
 **
 * Created: (28/11/2007) 
 *
 * Last Updated: (01/01/2008)   
 */

#ifndef _EXPMULTICLASSLOSS_CPP_
#define _EXPMULTICLASSLOSS_CPP_

#include <fstream>
#include <sstream>
#include "common.hpp"
#include "expmulticlassloss.hpp"
#include "configuration.hpp"

/**
 * Constructor
 *
 *   @param data [read/write] pointer to data container
 */
CExpMulticlassLoss::CExpMulticlassLoss(CModel* &model, CVecData* &data)
	: score(0),
          mygrad(0),
          val(0),
          f(0),
          l(0),
          mw(0),
          mgrad(0),
          _data(data)
{
        numOfW = 2;
	m = _data->slice_size();
	scalingFactor = 1.0/_data->size();
	
	// initialize w
	if(not model->IsInitialized())
	{                
		dimOfW = _data->dim();                
		model->Initialize(numOfW, dimOfW, _data->bias());
	}
	else
	{                
		dimOfW = model->GetDimOfW();
	}
	
        // MUST keep a pointer to the model
	_model = model;
	
        // determine the number of classes
        if(not data->HasLabel())
        {
                Configuration &config = Configuration::GetInstance();
                classes = config.GetInt("EXP_MULTICLASS.numOfClasses");
        }
        else
        {
                // sanity check
                assert((int)_data->MinLabel() > 0);
                classes = (size_t)_data->MaxLabel();                
        }
        
        // allocate memory for temporaries needed in loss and gradient computation
        score = new Scalar[classes];
	mygrad = new Scalar[numOfW];
        val = new Scalar[numOfW];        
        f = new TheMatrix(m, numOfW, SML::DENSE);        
        l = new TheMatrix(numOfW, m, SML::DENSE);
        mw = new TheMatrix(dimOfW, numOfW, SML::DENSE);
        mgrad = new TheMatrix(numOfW, dimOfW, SML::DENSE);
        
        if(verbosity)
        {
		std::cout << "In CExpMulticlassLoss::CExpMulticlassLoss() for multiclass" << std::endl;
                std::cout << "Number of classes: " << classes << std::endl;                
        }	
}

        
CExpMulticlassLoss::~CExpMulticlassLoss()
{
        if(score) delete[] score;
        if(mygrad) delete[] mygrad;        
        if(val) delete[] val;
        if(f) delete f;
        if(l) delete l;
        if(mw) delete mw;
        if(mgrad) delete mgrad;
}


/** 
 *  Display (usage) message. 
 */
 void CExpMulticlassLoss::Usage()
{
	std::cout << "Exponential Multiclass loss function:" << std::endl;
	std::cout << "No special parameters needed!" << std::endl;
}

/**  
 *  Compute loss function value. NOTE: We assume that the labels are
 *  integers which range from 1 to the data.maxLabel (typicall_y assumed
 *  to be small, of the order of 5 or so).
 *   
 *   @param loss [write] loss function value
 */
void CExpMulticlassLoss::ComputeLoss(Scalar& loss)
{
	TheMatrix &w = _model->GetW();	
	CalcF(w, *f);	
	Scalar* y = _data->labels().Data();	   
	loss = 0;
        
	for(unsigned int i = 0; i < m; i++)
	{
		unsigned int len = 0; 
		//Scalar val[numOfW];  //chteo: pushed into constructor
		f->GetRow(i, len, val);
		assert(len == (unsigned int)numOfW);
		CalcProb(val, score);
		// Loss is negative log-likelihood of the true label
		loss -= log(score[(unsigned int)y[i]-1]); 
	}
	
	//loss = loss/(1.0*m);  //chteo: should be divided by _data->size()
	loss = loss*scalingFactor;
}

/**  
 *  Compute loss function value and gradient of loss function w.r.t. w 
 *   
 *   @param loss [write] loss function value
 *   @param grad [write] gradient of loss function w.r.t w
 */
void CExpMulticlassLoss::ComputeLossAndGradient(Scalar& loss, TheMatrix& grad)
{
#if DEBUG
        std::cout << "numofw: " << numOfW << "   m: " << m << "  classes: " << classes << std::endl;
#endif        
	TheMatrix &w = _model->GetW();	
	CalcF(w, *f);	
	Scalar* y = _data->labels().Data();	
	loss = 0;
        
	for(unsigned int i = 0; i < m; i++)
	{
		unsigned int len = 0; 
                //Scalar val[numOfW];  //chteo: pushed into constructor
		f->GetRow(i, len, val);
		assert(len == (unsigned int)numOfW);
		CalcProb(val, score);
		
		// Gradient calculation
		mygrad[0] = -y[i]; 
		mygrad[1] = -(y[i]*y[i]);
		for(unsigned int k = 0; k < classes; k++)
		{
			mygrad[0] += score[k]*(k+1);
			mygrad[1] += score[k]*(k+1)*(k+1);
		}
		l->Set(0, i, mygrad[0]);
		l->Set(1, i, mygrad[1]);
                //chteo: gradient computation for multiple hyperplane can be sped up [see CMulticlassLoss for example]
                
		// Loss is negative log-likelihood of the true label
		loss -= log(score[(unsigned int)y[i]-1]); 
	}
	
	//loss = loss/(1.0*m);
	loss = loss*scalingFactor;
   
	// Matrix version of grad 	
	gradtime.Start();
	_data->Grad(*l, *mgrad);
	gradtime.Stop();
	
	// Flatten the rows of mgrad into a vector
	grad.Assign(*mgrad);
	
	//chteo: why there is no scaling for gradient vector? i add it here
	grad.Scale(scalingFactor);
}

void CExpMulticlassLoss::CalcProb(const Scalar* fval, Scalar* prob)
{
	// svnvish: BUGBUG 
	// In this function we make the explicit assumption that numOfW = 2 
	
	Scalar max_t = fval[0] + fval[1];
	for(unsigned int k = 0; k < classes; k++){
		prob[k] = fval[0]*(k+1) + fval[1]*(k+1)*(k+1);
		if(prob[k] > max_t)
			max_t = prob[k];
	}
	Scalar log_part = 0;
	for(unsigned int k = 0; k < classes; k++){
		prob[k] = exp(prob[k] - max_t);
		log_part += prob[k];
	}
	
	for(unsigned int k = 0; k < classes; k++){
		prob[k] /= log_part;
	}
}

void CExpMulticlassLoss::CalcF(const TheMatrix& w, TheMatrix& f)
{
	// Reshape w
	// The first column of mw contains w_{a} 
	// The second column of mw contains w_{b} 	
	
	// svnvish: BUGBUG
	// Explicitly assigning individual entries. Need to replace this.  		
        Scalar* w_array = w.Data();
	for(unsigned int i = 0; i < dimOfW; i++)
		mw->Set(i, 0, w_array[i]);
	
	for(unsigned int i = 0; i < dimOfW; i++)
		mw->Set(i, 1, w_array[i+dimOfW]);

	// compute f (size: m x 2). Since
	// <\phi(x, y), w> = <\phi(x), w_{a}> y + <\phi(x), w_{b}>y^2.
	// The first column of f contains entries <\phi(x), w_{a}>
	// The second column of f contains entries <\phi(x), w_{b}>
	_data->ComputeF(*mw, f);        
}


// chteo: will clean this up
void CExpMulticlassLoss::Predict(CModel* model)
{
        // sanity check
	assert(_data);
        
	TheMatrix &w = model->GetW();
        CalcF(w, *f);
        assert(f->Cols() == numOfW);
        
        // number of classes if a big problem in Predict() where dataset does not contain labels!
        // since we couldn't infer it from the number of W
        // so we read it from user's configuration file (this is done in the constructor)
                
        int* ybar = new int[m];
        Scalar* f_ybar = new Scalar[m];
        
        // argmax
	for(unsigned int i = 0; i < m; i++)
        {
		unsigned int len = 0; 
		f->GetRow(i, len, val);		
                
		CalcProb(val, score);		
                
		f_ybar[i] = score[0];
		ybar[i] = 1;
		for(unsigned int k = 1; k < classes; k++)
                {
			if(score[k] > f_ybar[i])
                        {
                                f_ybar[i] = score[k];
				ybar[i] = k+1;
                        }
		}
        }
        
        // write output to files
        Configuration &config = Configuration::GetInstance();
        if(config.GetString("Program.mode") != "LEARNING")
	{
                bool success = false;
                std::string predictedLabelsFn = "predictedLabels";
                if(config.IsSet("Prediction.predictedLabelsFile"))
                        predictedLabelsFn = config.GetString("Prediction.predictedLabelsFile");
                
                success = WriteFile<int*>(predictedLabelsFn, m, ybar);
                if(success) std::cout << "Predicted labels file written." << std::endl;
                
                std::string decisionFunctionValuesFn = "decisionFunctionValuesFile";
                if(config.IsSet("Prediction.decisionFunctionValuesFile"))
                        decisionFunctionValuesFn = config.GetString("Prediction.decisionFunctionValuesFile");
                
                success = WriteFile<Scalar*>(decisionFunctionValuesFn, m, f_ybar);
                if(success) std::cout << "Decision function values file written." << std::endl;
        }
        
	// clean up
	if(ybar) delete [] ybar;
	if(f_ybar) delete [] f_ybar;
}


void CExpMulticlassLoss::Evaluate(CModel* model)
{
      // sanity check
	assert(_data);
        
	TheMatrix &w = model->GetW();
        CalcF(w, *f);
        assert(f->Cols() == numOfW);
                
        int* ybar = new int[m];
        Scalar* f_ybar = new Scalar[m];
        
        // argmax
	for(unsigned int i = 0; i < m; i++)
        {
		unsigned int len = 0; 
		f->GetRow(i, len, val);		
                
		CalcProb(val, score);		
                
		f_ybar[i] = score[0];
		ybar[i] = 1;
		for(unsigned int k = 1; k < classes; k++)
                {
			if(score[k] > f_ybar[i])
                        {
                                f_ybar[i] = score[k];
				ybar[i] = k+1;
                        }
		}
        }
        
        // performance
        const TheMatrix& Y = _data->labels();
        int errorCnt = 0;
        for(unsigned int i=0; i<m; i++)
        {
                Scalar y_i = -1;
                Y.Get(i,y_i);
                if(ybar[i] != (int)y_i)
                        errorCnt++;
        }
        
        Scalar loss = 0.0;
        ComputeLoss(loss);
        
        std::cout << "\nPerformance:" << std::endl;
	std::cout << "1. Accuracy (%): " << 100.0*(m-errorCnt)/m << std::endl;
	std::cout << "2. Misclassification: " << errorCnt << std::endl;
        std::cout << "3. Loss: " << loss << std::endl;
        
        
        // write output to files
        Configuration &config = Configuration::GetInstance();
        if(config.GetString("Program.mode") != "LEARNING")
	{
                bool success = false;
                std::string predictedLabelsFn = "predictedLabels";
                if(config.IsSet("Prediction.predictedLabelsFile"))
                        predictedLabelsFn = config.GetString("Prediction.predictedLabelsFile");
                
                success = WriteFile<int*>(predictedLabelsFn, m, ybar);
                if(success) std::cout << "Predicted labels file written." << std::endl;
                
                std::string decisionFunctionValuesFn = "decisionFunctionValuesFile";
                if(config.IsSet("Prediction.decisionFunctionValuesFile"))
                        decisionFunctionValuesFn = config.GetString("Prediction.decisionFunctionValuesFile");
                
                success = WriteFile<Scalar*>(decisionFunctionValuesFn, m, f_ybar);
                if(success) std::cout << "Decision function values file written." << std::endl;
        }
        
	// clean up
	if(ybar) delete [] ybar;
	if(f_ybar) delete [] f_ybar;	
}

#endif
