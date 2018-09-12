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
 * Created: (26/01/2008) 
 *
 * Last Updated:
 */

#ifndef _MULTILABELLOSS_CPP_
#define _MULTILABELLOSS_CPP_

#include "common.hpp"  // WriteFile()
#include "multilabelloss.hpp"
#include "configuration.hpp"
#include <sstream>

/** Constructor
 *
 *  @param w    [r/w] external pointer to weight vector
 *  @param data [r/w] pointer to data container
 *  @param additiveLabelLoss [read] flag to determine whether to "add" label loss to loss function or not (default: true)
 */
CMultilabelLoss::CMultilabelLoss(CModel* &model, CMultilabelVecData* &data)
        : _data(data),
          f(0),
          matW(0),
          matG(0),
          g(0)          
{
	// initialize model
	if(not model->IsInitialized())
	{
		if(not data->HasLabel())
			throw CBMRMException("Error: dataset contains no labels; unable to determine number of classes!",
								 "CMultilabelLoss::CMultilabelLoss()");
		numOfClass = (int)data->MaxLabel();
		model->Initialize(numOfClass, _data->dim(), _data->bias());
	}
	else
	{
		numOfClass = model->GetNumOfW();
	}

	// keep a pointer to model
	_model = model;
		
        // check if labels are suitable for multilabel loss
	if(_data->HasLabel() && (_data->MinLabel() < 1.0)) 
        	throw CBMRMException("Labels must be positive integers","CMultilabelLoss::CMultilabelLoss()");
	
	// check if we can access the examples as individual vectors
	if(! _data->HasMatrixRowView()) 
        {
                if(verbosity)
                        std::cout << "CMultilabelLoss: creating explicit views to examples... " << std::endl;
                _data->CreateFeatureMatrixRowViews();
	}	
	                	
	// determine the number of classes from dataset
        m = _data->slice_size();
        scalingFactor = scalingFactor/_data->size();
        matW = new TheMatrix(numOfClass, _data->dim(), SML::DENSE);
        matG = new TheMatrix(numOfClass, _data->dim(), SML::DENSE);
        f = new TheMatrix(1, numOfClass, SML::DENSE);
		
        //i.e. g[i] is the i-th row of matG
        g = (TheMatrix**)malloc(sizeof(TheMatrix*)*numOfClass);
        for(unsigned int i=0; i<numOfClass; i++)
                g[i] = matG->CreateMatrixRowView(i);
      
	if(verbosity > 0)
	{
                std::cout << "In CMultilabelLoss::CMultilabelLoss()" << std::endl;
                std::cout << "Problem: " << numOfClass << "-class multilabel classification" << std::endl;
	}
}


/**  Destructor 
 */
CMultilabelLoss::~CMultilabelLoss()
{
	if(matW) {delete matW; matW = 0;}
	if(matG) {delete matG; matG = 0;}
	if(f) {delete f; f = 0;}
	if(g) 
        { 
            for(unsigned int i=0; i<numOfClass; i++)
                delete g[i];            
            free(g); 
            g = 0;
        }
}


/** 
 * Display (usage) message. 
 */
void CMultilabelLoss::Usage()
{
	std::cout << "Max-margin multilabel loss function:" << std::endl;
        std::cout << "  No parameter needed" << std::endl;
}


/**  Compute loss
 */
void CMultilabelLoss::ComputeLoss(Scalar& loss)
{
	TheMatrix &w = _model->GetW();
	//int min_y_i = 0;        
	double min_f_y_i = 0;
        //int max_ybar_i = 0;
	double max_f_ybar_i = 0;
	
	// assign the content of w (1D) to matW (matrix)
	matW->Assign(w);
	loss = 0;
	const vector<vector<Scalar> > &Y = _data->labels();
	
	for(unsigned int i=0; i<m; i++) 
	{
                //min_y_i = 0;
                min_f_y_i = SML::INFTY;
                //max_ybar_i = 0;
                max_f_ybar_i = -SML::INFTY;
                
		_data->ComputeFi(*matW, *f, i);
                // look for min_{y} f(x,y) and max_{y'} f(x,y'), where y is in true label set and y' in the complement
                // loss := max(0, 1 + max_{y'} f(x,y') - min_{y} f(x,y)) 
                unsigned int idx = 0;
                Scalar f_k = 0;
                for(unsigned int k=0; k<numOfClass; k++)
                {
                        f->Get(k, f_k);
                        
                        if((idx < Y[i].size()) && (k+1 == Y[i][idx]))
                        {
                                if(f_k <= min_f_y_i)
                                {
                                        min_f_y_i = f_k;
                                        //min_y_i = k+1;
                                }
                                idx++;                                        
                        }
                        else
                        {
                                if(f_k >= max_f_ybar_i)
                                {
                                        max_f_ybar_i = f_k;
                                        //max_ybar_i = k+1;
                                }
                        }
                }
                loss += (Scalar)std::max(0.0, 1.0 + max_f_ybar_i - min_f_y_i);
	}
	loss *= scalingFactor;  
}




/**   Compute loss and gradient
 *    Note: class label starts from 1
 */
void CMultilabelLoss::ComputeLossAndGradient(Scalar& loss, TheMatrix& grad)
{
	TheMatrix &w = _model->GetW();
	unsigned int min_y_i = 0;        
	double min_f_y_i = 0;
        unsigned int max_ybar_i = 0;
	double max_f_ybar_i = 0;
	
	// assign the content of w (1D) to matW (matrix)
	matW->Assign(w);  
	matG->Zero();
	loss = 0;
        const vector<vector<Scalar> > &Y = _data->labels();	
	
        for(unsigned int i=0; i<m; i++) 
	{
                //printf("[%d]",i); fflush(stdout);
                min_y_i = 1;
                min_f_y_i = SML::INFTY;
                max_ybar_i = 1;
                max_f_ybar_i = -SML::INFTY;
                
		_data->ComputeFi(*matW, *f, i);
                // look for min_{y} f(x,y) and max_{y'} f(x,y'), where y is in true label set and y' in the complement
                // loss := max(0, 1 + max_{y'} f(x,y') - min_{y} f(x,y)) 
                unsigned int idx = 0;
                Scalar f_k = 0;
               for(unsigned int k=0; k<numOfClass; k++)
               {
                       //printf("."); fflush(stdout);
                       f->Get(k, f_k);
                       
                       if((idx < Y[i].size()) && (k+1 == Y[i][idx]))
                       {
                               if(f_k <= min_f_y_i)
                               {
                                       min_f_y_i = f_k;
                                       min_y_i = k+1;
                               }
                               idx++;                                        
                       }
                       else
                       {
                               if(f_k >= max_f_ybar_i)
                               {
                                       max_f_ybar_i = f_k;
                                       max_ybar_i = k+1;
                               }
                       }
               }
//                 for(unsigned int k=0; k<numOfClass; k++)
//                 {
//                         //printf("."); fflush(stdout);
//                         f->Get(k, f_k);
                        
//                         for(idx=0; idx<Y[i].size(); idx++)
//                         {
//                                 if(Y[i][idx] == k+1)
//                                 {
//                                         if(f_k <= min_f_y_i)
//                                         {
//                                                 min_f_y_i = f_k;
//                                                 min_y_i = k+1;
//                                         }                                
//                                 }
//                                 else
//                                 {
//                                         if(f_k >= max_f_ybar_i)
//                                         {
//                                                 max_f_ybar_i = f_k;
//                                                 max_ybar_i = k+1;
//                                         }
//                                 }

//                         }
//                 }

                // loss for example i
                loss += (Scalar)std::max(0.0, 1.0 + max_f_ybar_i - min_f_y_i);
                
                // gradient for example i
                _data->AddElement(*g[min_y_i-1], i, -1.0);                
                _data->AddElement(*g[max_ybar_i-1], i, 1.0);
	}
	
	grad.Assign(*matG);
	grad.Scale(scalingFactor);	
	loss *= scalingFactor;        
        
        //Evaluate(_model);
}



void CMultilabelLoss::Predict(CModel* model)
{
        // sanity check
        assert(model);
        assert(_data);
        assert(m > 0);
               
        // more defensive precaution
        assert(_data->dim() == model->GetDimOfW());
        
        if(numOfClass != model->GetNumOfW())
        {
                if(matW) 
                {
                        delete matW;                        
                        matW = new TheMatrix(model->GetNumOfW(), model->GetDimOfW(), SML::DENSE);
                }                        
        }
        
        numOfClass = std::min(numOfClass, model->GetNumOfW());                      
        assert(matW);
        
        int* ybar = new int[m];
	Scalar* f_ybar = new Scalar[m];
        
        DecisionAndPrediction(model, ybar, f_ybar);
 
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
 


/**  Evaluate the performance of the model on the examples in _data
 *
 *   @param model [read] model object containing a weight vector
 */

void CMultilabelLoss::Evaluate(CModel* model)
{
        // sanity check
        assert(model);
        assert(_data);
        assert(m > 0);
        
        if(not _data->HasLabel())
        {
                throw CBMRMException("Data object does not contain labels","CMultilabelLoss::Evaluate()");
        }
        
        // more defensive precaution
        assert(_data->HasLabel());
        assert(_data->dim() == model->GetDimOfW());
        
        if(numOfClass != model->GetNumOfW())
        {
                if(matW) 
                {
                        delete matW;                        
                        matW = new TheMatrix(model->GetNumOfW(), model->GetDimOfW(), SML::DENSE);
                }                        
        }
        
        numOfClass = std::min(numOfClass, model->GetNumOfW());                      
        assert(matW);
        
        int* ybar = new int[m];
	Scalar* f_ybar = new Scalar[m];
        
        DecisionAndPrediction(model, ybar, f_ybar);

        // evaluation	
        int oneErrorCnt = m;
        const vector<vector<Scalar> > &Y = _data->labels();
        
        for(unsigned int i=0; i < m; i++)
        {
                for(unsigned int j=0; j<Y[i].size(); j++)
                        if(ybar[i] == Y[i][j])
                        {
                                oneErrorCnt--;
                                break;
                        }
        }
        
        Scalar loss = 0.0;
        ComputeLoss(loss);
		
        // performance
	std::cout << "\nPerformance:" << std::endl;
	std::cout << "1. OneErr (%): " << 100.0*(oneErrorCnt)/m << std::endl;
        std::cout << "2. Acc@1 (%): " << 100.0*(m-oneErrorCnt)/m << std::endl;
        std::cout << "3. Loss: " << loss << std::endl;
                
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


void CMultilabelLoss::DecisionAndPrediction(CModel* model, int *ybar, Scalar *f_ybar)
{
        // sanity check	
        assert(ybar);
        assert(f_ybar);
        
        TheMatrix &w = model->GetW();        
        matW->Assign(w);        	
	        
	for(unsigned int i=0; i < m; i++)
	{
		_data->ComputeFi(*matW, *f, i);
		Scalar max_f_ybar_i = -SML::INFTY;
		
		for(unsigned int j = 0; j < numOfClass; j++)
		{
			f->Get(j, f_ybar[i]);
			if(f_ybar[i] > max_f_ybar_i)
			{
				max_f_ybar_i = f_ybar[i];
				ybar[i] = j+1;
			}
		} 		
	}                        
}

#endif
