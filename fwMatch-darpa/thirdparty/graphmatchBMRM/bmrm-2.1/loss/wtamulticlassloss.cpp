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

#ifndef _WTAMULTICLASSLOSS_CPP_
#define _WTAMULTICLASSLOSS_CPP_

#include "common.hpp"  // WriteFile()
#include "wtamulticlassloss.hpp"
#include "configuration.hpp"
#include <sstream>

/** Constructor
 *
 *  @param w    [r/w] external pointer to weight vector
 *  @param data [r/w] pointer to data container
 *  @param additiveLabelLoss [read] flag to determine whether to "add" label loss to loss function or not (default: true)
 */
CWTAMulticlassLoss::CWTAMulticlassLoss(CModel* &model, CVecData* &data, bool additiveLabelLoss)
        : _data(data),
          f(0),
          matW(0),
          matG(0),
          g(0),
          _additiveLabelLoss(additiveLabelLoss)
          
{
	// initialize model
	if(not model->IsInitialized())
	{
		if(not data->HasLabel())
			throw CBMRMException("Error: dataset contains no labels; unable to determine number of classes!",
								 "CWTAMulticlassLoss::CWTAMulticlassLoss()");
		numOfClass = (int)data->MaxLabel();
		model->Initialize(numOfClass, _data->dim(), _data->bias());
	}
	else
	{
		numOfClass = model->GetNumOfW();
	}

	// keep a pointer to model
	_model = model;
	
	
        // check if labels are suitable for multiclass classification
	if(_data->HasLabel() and ((int)_data->MinLabel() < 1.0)) {
		throw CBMRMException("Labels must be positive integers",
							 "CWTAMulticlassLoss::CWTAMulticlassLoss()");
	}
	
	// check if we can access the examples as individual vectors
	if(not _data->HasMatrixRowView()) {
                if(verbosity)
                        std::cout << "CWTAMulticlassLoss: creating explicit views to examples... " << std::endl;
                _data->CreateFeatureMatrixRowViews();
                _data->CreateLabelMatrixRowViews();		
	}
	
	
	Configuration &config = Configuration::GetInstance();
	
        // default margin scaling function
        Gamma = &CWTAMulticlassLoss::MarginRescaling;
        marginScalingFunctionType = "MARGIN_RESCALING";
        
	if(config.IsSet("WTA_MULTICLASS.marginScalingFunctionType"))
        {
		if(config.GetString("WTA_MULTICLASS.marginScalingFunctionType") == "SLACK_RESCALING")
                {
                        Gamma = &CWTAMulticlassLoss::SlackRescaling;
                        marginScalingFunctionType = "SLACK_RESCALING";
                }
        }
        
	// default label loss function
        Delta = &CWTAMulticlassLoss::ZeroOne;
        labelLossType = "ZERO_ONE";

        if(config.IsSet("WTA_MULTICLASS.labelLossType"))
        {
                if(config.GetString("WTA_MULTICLASS.labelLossType") == "SQUARED_DIFFERENCE")
                {
                        Delta = &CWTAMulticlassLoss::SquaredDifference;
                        labelLossType = "SQUARED_DIFFERENCE";
                }
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
                std::cout << "In CWTAMulticlassLoss::CWTAMulticlassLoss()" << std::endl;
                std::cout << "Problem: " << numOfClass << "-class classification" << std::endl;
                std::cout << "Margin scaling function type: " << marginScalingFunctionType << std::endl;
                std::cout << "Label loss type: " << labelLossType << std::endl << std::endl;                
	}
}


/**  Destructor 
 */
CWTAMulticlassLoss::~CWTAMulticlassLoss()
{
	if(matW) {delete matW; matW = 0;}
	if(matG) {delete matG; matG = 0;}
	if(f) {delete f; f = 0;}
	if(g) { 
            for(unsigned int i=0; i<numOfClass; i++) {
                delete g[i];
            }
            free(g); 
            g = 0;
        }
}


/** 
 * Display (usage) message. 
 */
void CWTAMulticlassLoss::Usage()
{
	std::cout << "Winner-takes-all multiclass loss function:" << std::endl;
        std::cout << "1. Types of loss rescaling (String Multiclass.rescaling):" << std::endl;
        std::cout << "* Slack re-scaling (SLACK_RESCALING) -- default" << std::endl;
        std::cout << "* Margin re-scaling (MARGIN_RESCALING)" << std::endl;
        std::cout << "2. Label loss type (string Multiclass.labelLossType):  " << std::endl;
	std::cout << "* 0/1 (ZEROONE) -- default" << std::endl;
	std::cout << "* Squared difference (SQUARED_DIFFERENCE)" << std::endl;
}


/**  Compute loss
 */
void CWTAMulticlassLoss::ComputeLoss(Scalar& loss)
{
	TheMatrix &w = _model->GetW();
	Scalar y_i = 1;
	Scalar ybar_i = 1;
	Scalar f_y_i = 0;
	Scalar f_ybar_i = 0;
	Scalar loss_i = 0;
	Scalar max_loss_i = 0;
	Scalar labelloss = 0;
	Scalar max_labelloss = 0;
	
	// assign the content of w (1D) to matW (matrix)
	matW->Assign(w);  // chteo: to be replace with TheMatrix::StackRow()
	loss = 0;
	const TheMatrix &Y = _data->labels();
	
	for(unsigned int i=0; i<m; i++) 
	{
		_data->ComputeFi(*matW, *f, i);
		Y.Get(i, y_i);
		f->Get((int)y_i-1, f_y_i);
		max_loss_i = 0.0;
		loss_i = 0.0;
		ybar_i = -1;
		
		for(unsigned int j=0; j<numOfClass; j++) 
                {
			f->Get(j, f_ybar_i);
			labelloss = (this->*Delta)((int)y_i, j+1);		
                        loss_i = (this->*Gamma)(labelloss) * (f_ybar_i - f_y_i); // + labelloss
                        
                        // this is modified for ramp loss!
                        if(_additiveLabelLoss)
                                loss_i += labelloss;
                        			
			if(loss_i > max_loss_i) {
				max_loss_i = loss_i;
				max_labelloss = labelloss;
				ybar_i = j+1;
			}
		}
		
		if(ybar_i > 0 and ((int)ybar_i != (int)y_i)) 
                {
			loss += max_loss_i;
		}
	}
	loss *= scalingFactor;  
}




/**   Compute loss and gradient
 *    Note: class label starts from 1
 */
void CWTAMulticlassLoss::ComputeLossAndGradient(Scalar& loss, TheMatrix& grad)
{
	TheMatrix &w = _model->GetW();
	Scalar y_i = 1;
	Scalar ybar_i = 1;
	Scalar f_y_i = 0;
	Scalar f_ybar_i = 0;
	Scalar loss_i = 0;
	Scalar max_loss_i = 0;
	Scalar labelloss = 0;
        Scalar marginscale = 0;
	Scalar max_marginscale = 0;
	
	// assign the content of w (1D) to matW (matrix)
	matW->Assign(w);  
	matG->Zero();
	loss = 0;
	const TheMatrix &Y = _data->labels();
	
        int errcnt = 0;
	argmaxtime.Start();
	for(unsigned int i=0; i<m; i++) 
	{
                _data->ComputeFi(*matW, *f, i);
                Y.Get(i, y_i);
                f->Get((int)y_i-1, f_y_i);
                max_loss_i = 0.0;
                loss_i = 0.0;
                ybar_i = -1;
		
                // Do: argmax_y Gamma(Delta(y,ybar))(\f_ybar - f_y) + Delta(y,ybar), where f_* := <w,*>
                for(unsigned int j=0; j<numOfClass; j++) 
                {
                        f->Get(j, f_ybar_i);
                        labelloss = (this->*Delta)((int)y_i, j+1);
			marginscale = (this->*Gamma)(labelloss);                               		
                        loss_i = marginscale * (f_ybar_i - f_y_i); // + labelloss;
                      
                        // this is modified for ramp loss
                        if(_additiveLabelLoss)
                                loss_i += labelloss;
                                
                        if(loss_i > max_loss_i) 
                        {
                                max_loss_i = loss_i;
                                max_marginscale = marginscale;
                                ybar_i = j+1;
                        }
                }
                gradtime.Start();
                if(ybar_i > 0 and ((int)ybar_i != (int)y_i)) 
                {
                        loss += max_loss_i;
                        _data->AddElement(*g[(int)y_i-1], i, -max_marginscale);
                        _data->AddElement(*g[(int)ybar_i-1], i, max_marginscale);
                        errcnt++;
                }
                gradtime.Stop();
        }
	argmaxtime.Stop();
	
	gradtime.Start();
	grad.Assign(*matG);
	grad.Scale(scalingFactor);
        gradtime.Stop();
	
	loss *= scalingFactor;
        
	if(verbosity > 0) 
	{
                std::cout << std::endl;
                std::cout << "argmax + w*x + fill l time: " << argmaxtime.CPUTotal() - gradtime.CPUTotal() << std::endl;
                std::cout << "X*l time: " << gradtime.CPUTotal() << std::endl;
                std::cout << "errcnt: " << errcnt << std::endl;
	}
}




void CWTAMulticlassLoss::Predict(CModel* model)
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

void CWTAMulticlassLoss::Evaluate(CModel* model)
{
        // sanity check
        assert(model);
        assert(_data);
        assert(m > 0);
        
        if(not _data->HasLabel())
        {
                throw CBMRMException("Data object does not contain labels","CWTAMulticlassLoss::Evaluate()");
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
        int errorCnt = 0;
        Scalar y_i = 0;;
        const TheMatrix &Y = _data->labels();
        
        for(unsigned int i=0; i < m; i++)
        {
                Y.Get(i, y_i);
                if(ybar[i] != (int)y_i) 
			errorCnt++;
        }
        
        Scalar loss = 0.0;
        ComputeLoss(loss);
		
        // performance
	std::cout << "\nPerformance:" << std::endl;
	std::cout << "1. Accuracy (%): " << 100.0*(m-errorCnt)/m << std::endl;
	std::cout << "2. Misclassification: " << errorCnt << std::endl;
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


void CWTAMulticlassLoss::DecisionAndPrediction(CModel* model, int *ybar, Scalar *f_ybar)
{
        // sanity check	
        assert(ybar);
        assert(f_ybar);
        
        TheMatrix &w = model->GetW();        
        matW->Assign(w);        	
	        
	for(unsigned int i=0; i < m; i++)
	{
	  //Scalar y_i = 0.0;		
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
